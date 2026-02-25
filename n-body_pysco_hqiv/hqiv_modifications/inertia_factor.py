"""
HQIV Inertia Factor Module
==========================

Computes the thermodynamic inertia reduction factor f(a_loc, φ) from the action principle.

Paper Reference: paper/main.tex, Section 4.2
    f(a_loc, φ) = max( a_loc / (a_loc + cφ/6), f_min )

where:
    - a_loc is the local acceleration magnitude
    - φ is the horizon field
    - cφ/6 gives the minimum acceleration scale
    - f_min ≈ 0.01 is the thermodynamic floor

The factor χ ≈ 0.172 from the full light-cone average (Brodie's overlap integral)
relates the minimum acceleration to the horizon field: a_min = χ c φ

Key Physics:
------------
- High acceleration (a_loc >> a_min): f → 1 (standard inertia)
- Low acceleration (a_loc ~ a_min): f < 1 (reduced inertia)
- The floor f_min prevents breakdown in deep MOND regime

Author: HQIV Team
"""

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
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


# Load default parameters from config
_config = load_config()
_chi = float(_config.get('hqiv_parameters', 'chi', fallback='0.172'))
_f_min = float(_config.get('hqiv_parameters', 'f_min', fallback='0.01'))
_inertia_form = _config.get('hqiv_parameters', 'inertia_form', fallback='thermo')


def inertia_reduction_factor(alpha, phi, chi=_chi, f_min=_f_min, form=_inertia_form):
    """
    Compute the inertia reduction factor f(α, φ) with thermodynamic floor.
    
    Paper Reference: paper/main.tex, Eq. in Section 4.2
        f(α, φ) = max( α / (α + c φ / 6), f_min )
    
    This is the thermodynamic form from Brodie's overlap integral.
    
    Parameters
    ----------
    alpha : float or array
        Local acceleration magnitude |a_local| [m/s²]
    phi : float or array
        Local horizon field φ [m/s²] (same units as alpha)
    chi : float
        Scaling factor from light-cone average (default 0.172)
    f_min : float
        Thermodynamic floor (default 0.01)
    form : str
        Form of inertia reduction:
        - 'thermo': Thermodynamic form (recommended, from Brodie)
        - 'sqrt': Square-root form from earlier QI
        - 'brodie': Brodie's linear form without floor
        
    Returns
    -------
    f : float or array
        Inertia reduction factor (f_min ≤ f ≤ 1)
        
    Notes
    -----
    The thermodynamic form is derived from:
        f(a, Θ) = a / (a + cH/6)
    
    where the factor 1/6 comes from the backward-hemisphere overlap integral:
        ∫₀^{π/2} cos²θ sinθ dθ = 1/6
    
    With χ = 0.172 from the full light-cone average, the minimum acceleration is:
        a_min = χ c φ
    
    Examples
    --------
    >>> # High acceleration regime (standard inertia)
    >>> f = inertia_reduction_factor(1e-9, 1e-10)  # a >> a_min
    >>> print(f)  # Should be close to 1
    
    >>> # Low acceleration regime (reduced inertia)
    >>> f = inertia_reduction_factor(1e-12, 1e-10)  # a ~ a_min
    >>> print(f)  # Should be < 1
    """
    alpha = np.asarray(alpha)
    phi = np.asarray(phi)
    
    if form == 'thermo':
        # Thermodynamic floor from Brodie's thermodynamics
        # a_min = χ c φ (minimum acceleration from horizon information)
        # The factor cφ/6 comes from the overlap integral
        a_min = chi * c * phi
        
        # Brodie's linear form with thermodynamic floor
        # f = α / (α + a_min/6) = α / (α + χ c φ / 6)
        denominator = alpha + a_min / 6.0
        
        # Avoid division by zero
        denominator = np.maximum(denominator, 1e-100)
        
        f = alpha / denominator
        
        # Apply thermodynamic floor
        f = np.maximum(f, f_min)
        
        # Also cap at 1 (shouldn't exceed standard inertia)
        f = np.minimum(f, 1.0)
        
        return f
    
    elif form == 'sqrt':
        # HQIV square-root form (paper after eq. 8)
        # f = √(1 - cφ/α)
        ratio = c * phi / (alpha + 1e-100)
        
        # Clamp to avoid negative values under sqrt
        ratio = np.minimum(ratio, 0.99)
        
        f = np.sqrt(1.0 - ratio)
        
        # Apply floor
        f = np.maximum(f, f_min)
        
        return f
    
    elif form == 'brodie':
        # Brodie's form without floor
        f = alpha / (alpha + c * phi / 6.0)
        
        # Still apply floor for numerical stability
        f = np.maximum(f, f_min)
        
        return f
    
    else:
        raise ValueError(f"Unknown form: {form}. Use 'thermo', 'sqrt', or 'brodie'")


def compute_local_acceleration(acceleration_field, positions, dx, N_grid):
    """
    Compute the local acceleration magnitude at particle positions.
    
    Parameters
    ----------
    acceleration_field : array
        Gravitational acceleration field [3, N, N, N] in code units
    positions : array
        Particle positions [N_part, 3] in comoving coordinates
    dx : float
        Grid spacing
    N_grid : int
        Grid size
        
    Returns
    -------
    alpha : array
        Local acceleration magnitude at each particle position [code units]
    """
    # Interpolate acceleration to particle positions using CIC
    acc_p = np.zeros((len(positions), 3))
    
    for i, pos in enumerate(positions):
        ix = int(np.floor(pos[0] / dx)) % N_grid
        iy = int(np.floor(pos[1] / dx)) % N_grid
        iz = int(np.floor(pos[2] / dx)) % N_grid
        
        ux = pos[0] / dx - np.floor(pos[0] / dx)
        uy = pos[1] / dx - np.floor(pos[1] / dx)
        uz = pos[2] / dx - np.floor(pos[2] / dx)
        
        for ddx in [0, 1]:
            for ddy in [0, 1]:
                for ddz in [0, 1]:
                    w = ((1-ux) if ddx==0 else ux) * \
                        ((1-uy) if ddy==0 else uy) * \
                        ((1-uz) if ddz==0 else uz)
                    acc_p[i, 0] += w * acceleration_field[0, (ix+ddx)%N_grid, (iy+ddy)%N_grid, (iz+ddz)%N_grid]
                    acc_p[i, 1] += w * acceleration_field[1, (ix+ddx)%N_grid, (iy+ddy)%N_grid, (iz+ddz)%N_grid]
                    acc_p[i, 2] += w * acceleration_field[2, (ix+ddx)%N_grid, (iy+ddy)%N_grid, (iz+ddz)%N_grid]
    
    # Compute magnitude
    alpha = np.sqrt(np.sum(acc_p**2, axis=1))
    
    return alpha


def directional_inertia_factor(alpha, phi, grad_phi, velocity, 
                               chi=_chi, f_min=_f_min, form=_inertia_form):
    """
    Compute direction-dependent inertia factor.
    
    Paper Reference: paper/main.tex, discussion of directional dependence
        The inertia modification depends on the alignment between the local
        acceleration and the horizon gradient direction.
    
    For gas particles (high thermal acceleration): f → 1
    For galaxy particles (low acceleration): f < 1 with directional dependence
    
    Parameters
    ----------
    alpha : float or array
        Local acceleration magnitude
    phi : float or array
        Horizon field φ
    grad_phi : array
        Gradient of φ [3, ...] or [3] for single point
    velocity : array
        Particle velocity [3, ...] or [3] for single point
    chi : float
        Scaling factor
    f_min : float
        Thermodynamic floor
    form : str
        Form of inertia reduction
        
    Returns
    -------
    f : float or array
        Direction-dependent inertia factor
    f_parallel : float or array
        Factor for motion parallel to ∇φ
    f_perp : float or array
        Factor for motion perpendicular to ∇φ
        
    Notes
    -----
    The directional dependence arises because the horizon modification
    affects motion differently along vs. perpendicular to the horizon gradient.
    
    For the Bullet Cluster, this creates the observed separation between
    gas (collisional, high acceleration, f≈1) and galaxies (collisionless,
    low acceleration, f<1).
    """
    # Base inertia factor
    f_base = inertia_reduction_factor(alpha, phi, chi, f_min, form)
    
    # Compute direction of horizon gradient
    grad_phi_mag = np.sqrt(np.sum(grad_phi**2, axis=0)) + 1e-100
    
    # Unit vector along ∇φ
    e_grad = grad_phi / grad_phi_mag
    
    # Project velocity onto gradient direction
    v_parallel = np.sum(velocity * e_grad, axis=0)
    
    # Directional enhancement factor
    # Motion along ∇φ experiences different inertia than perpendicular
    # This is a simplified model; full treatment requires the action
    v_mag = np.sqrt(np.sum(velocity**2, axis=0)) + 1e-100
    
    # Cosine of angle between velocity and ∇φ
    cos_theta = v_parallel / v_mag
    cos_theta = np.clip(cos_theta, -1, 1)
    
    # Directional modification (simplified model)
    # Parallel motion: slightly enhanced inertia reduction
    # Perpendicular motion: standard reduction
    directional_factor = 1.0 + 0.1 * cos_theta**2
    
    # Apply directional factor
    f = f_base * directional_factor
    f = np.clip(f, f_min, 1.0)
    
    # Return components for analysis
    f_parallel = f_base * (1.0 + 0.1)  # Enhanced along gradient
    f_perp = f_base * (1.0 - 0.05)     # Slightly reduced perpendicular
    
    f_parallel = np.clip(f_parallel, f_min, 1.0)
    f_perp = np.clip(f_perp, f_min, 1.0)
    
    return f, f_parallel, f_perp


def inertia_factor_for_species(alpha, phi, species_type, 
                                chi=_chi, f_min=_f_min, form=_inertia_form):
    """
    Compute inertia factor for different particle species.
    
    For the Bullet Cluster simulation, we have:
    - Gas particles: High thermal acceleration → f ≈ 1
    - Galaxy particles: Low acceleration → f < 1
    
    Parameters
    ----------
    alpha : float or array
        Local acceleration magnitude
    phi : float or array
        Horizon field
    species_type : str or array
        Particle type: 'gas' or 'galaxy'
    chi : float
        Scaling factor
    f_min : float
        Thermodynamic floor
    form : str
        Form of inertia reduction
        
    Returns
    -------
    f : float or array
        Inertia factor for each particle
        
    Notes
    -----
    Gas particles have high thermal accelerations from the hot intracluster
    medium (T ~ 10^8 K), so they naturally have f → 1.
    
    Galaxy particles are collisionless and have low internal accelerations,
    so they experience the full inertia reduction.
    """
    # Base inertia factor
    f_base = inertia_reduction_factor(alpha, phi, chi, f_min, form)
    
    if isinstance(species_type, str):
        # Single species type
        if species_type == 'gas':
            # Gas has high thermal acceleration, f should be close to 1
            # Add thermal acceleration contribution
            # For ICM at T ~ 10^8 K, thermal acceleration ~ kT / (m_p * r)
            # This is typically >> a_min
            thermal_boost = 10.0  # Factor to boost effective acceleration
            f = inertia_reduction_factor(alpha * thermal_boost, phi, chi, f_min, form)
        else:
            f = f_base
    else:
        # Array of species types
        f = np.where(
            species_type == 'gas',
            inertia_reduction_factor(alpha * 10.0, phi, chi, f_min, form),
            f_base
        )
    
    return f


class InertiaFactorSolver:
    """
    Class to manage inertia factor computation throughout the simulation.
    
    Parameters
    ----------
    H0 : float
        Hubble constant [km/s/Mpc]
    chi : float
        Scaling factor from light-cone average
    f_min : float
        Thermodynamic floor
    form : str
        Form of inertia reduction
    """
    
    def __init__(self, H0=73.2, chi=_chi, f_min=_f_min, form=_inertia_form):
        self.H0 = H0 * 1000.0 / Mpc_m  # Convert to 1/s
        self.chi = chi
        self.f_min = f_min
        self.form = form
        
    def compute(self, alpha, phi):
        """Compute inertia factor for given acceleration and φ."""
        return inertia_reduction_factor(alpha, phi, self.chi, self.f_min, self.form)
    
    def compute_for_particles(self, acceleration, phi_particles, species=None):
        """
        Compute inertia factor for all particles.
        
        Parameters
        ----------
        acceleration : array
            Acceleration magnitude at each particle position
        phi_particles : array
            φ values at each particle position
        species : array, optional
            Species type for each particle ('gas' or 'galaxy')
            
        Returns
        -------
        f : array
            Inertia factor for each particle
        """
        if species is not None:
            return inertia_factor_for_species(
                acceleration, phi_particles, species,
                self.chi, self.f_min, self.form
            )
        else:
            return self.compute(acceleration, phi_particles)
    
    def apply_to_acceleration(self, acc, f):
        """
        Apply inertia factor to acceleration.
        
        Modified Euler equation (paper eq. 10):
            a_eff = a_grav / f
        
        Lower f means less inertia, so particles respond more strongly
        to gravitational acceleration.
        
        Parameters
        ----------
        acc : array
            Gravitational acceleration [N_part, 3]
        f : array
            Inertia factor [N_part]
            
        Returns
        -------
        acc_modified : array
            Modified acceleration [N_part, 3]
        """
        # Ensure f is properly shaped for broadcasting
        f = np.asarray(f).reshape(-1, 1)
        
        # Apply modification
        acc_modified = acc / f
        
        return acc_modified


# =============================================================================
# Numba-optimized functions for performance
# =============================================================================

try:
    from numba import njit, prange
    
    @njit(parallel=True, fastmath=True)
    def inertia_reduction_factor_numba(alpha, phi, chi, f_min):
        """
        Numba-optimized inertia factor computation.
        
        Uses the thermodynamic form for performance.
        """
        n = len(alpha)
        f = np.empty(n, dtype=np.float64)
        
        for i in prange(n):
            a_min = chi * c * phi[i]
            denom = alpha[i] + a_min / 6.0
            if denom > 0:
                f[i] = alpha[i] / denom
            else:
                f[i] = f_min
            
            # Apply floor
            if f[i] < f_min:
                f[i] = f_min
            if f[i] > 1.0:
                f[i] = 1.0
        
        return f
    
    _HAS_NUMBA = True
    
except ImportError:
    _HAS_NUMBA = False
    inertia_reduction_factor_numba = None


def compute_inertia_factor_fast(alpha, phi, chi=_chi, f_min=_f_min):
    """
    Fast inertia factor computation using Numba if available.
    
    Parameters
    ----------
    alpha : array
        Local acceleration magnitude
    phi : array
        Horizon field
    chi : float
        Scaling factor
    f_min : float
        Thermodynamic floor
        
    Returns
    -------
    f : array
        Inertia factor
    """
    if _HAS_NUMBA:
        return inertia_reduction_factor_numba(
            np.asarray(alpha, dtype=np.float64),
            np.asarray(phi, dtype=np.float64),
            chi, f_min
        )
    else:
        return inertia_reduction_factor(alpha, phi, chi, f_min, form='thermo')


# =============================================================================
# Testing
# =============================================================================

def test_inertia_factor():
    """Test the inertia factor computation."""
    import matplotlib.pyplot as plt
    
    print("Testing inertia factor computation...")
    
    # Test range of accelerations
    alpha = np.logspace(-14, -6, 1000)  # m/s²
    
    # Typical φ value (c * H0)
    H0 = 73.2 * 1000.0 / Mpc_m  # 1/s
    phi = c * H0  # m/s²
    
    print(f"φ = cH0 = {phi:.2e} m/s²")
    print(f"a_min = χ c φ = {0.172 * phi:.2e} m/s²")
    
    # Compute inertia factor for different forms
    f_thermo = inertia_reduction_factor(alpha, phi, form='thermo')
    f_sqrt = inertia_reduction_factor(alpha, phi, form='sqrt')
    f_brodie = inertia_reduction_factor(alpha, phi, form='brodie')
    
    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot f vs acceleration
    ax = axes[0]
    ax.loglog(alpha, f_thermo, 'b-', label='Thermodynamic', linewidth=2)
    ax.loglog(alpha, f_sqrt, 'r--', label='Square-root', linewidth=2)
    ax.loglog(alpha, f_brodie, 'g:', label='Brodie (no floor)', linewidth=2)
    
    # Mark characteristic scales
    a_min = 0.172 * phi
    ax.axvline(a_min, color='k', linestyle=':', alpha=0.5, label=f'a_min = {a_min:.2e}')
    ax.axhline(0.01, color='gray', linestyle='--', alpha=0.5, label='Floor f_min = 0.01')
    
    ax.set_xlabel('Local acceleration α [m/s²]')
    ax.set_ylabel('Inertia factor f')
    ax.set_title('HQIV Inertia Reduction Factor')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0.005, 1.5])
    
    # Plot f vs acceleration (linear scale for transition region)
    ax = axes[1]
    ax.semilogx(alpha, f_thermo, 'b-', label='Thermodynamic', linewidth=2)
    ax.semilogx(alpha, f_sqrt, 'r--', label='Square-root', linewidth=2)
    
    ax.axvline(a_min, color='k', linestyle=':', alpha=0.5)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3)
    
    ax.set_xlabel('Local acceleration α [m/s²]')
    ax.set_ylabel('Inertia factor f')
    ax.set_title('Inertia Factor Transition')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.1])
    
    plt.tight_layout()
    plt.savefig('inertia_factor_test.png', dpi=150)
    plt.close()
    print("Saved inertia_factor_test.png")
    
    # Print summary
    print("\nInertia Factor Summary:")
    print(f"  At a = a_min: f = {inertia_reduction_factor(a_min, phi):.3f}")
    print(f"  At a = 10*a_min: f = {inertia_reduction_factor(10*a_min, phi):.3f}")
    print(f"  At a = 0.1*a_min: f = {inertia_reduction_factor(0.1*a_min, phi):.3f}")
    
    return f_thermo


if __name__ == "__main__":
    test_inertia_factor()