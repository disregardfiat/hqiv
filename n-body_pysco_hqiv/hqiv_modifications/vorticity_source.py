"""
HQIV Vorticity Source Module
============================

Implements the action-derived vorticity source term from the HQIV framework.

Paper Reference: paper/main.tex, Section 6
    ∂ω/∂t + (v·∇)ω = (∂f/∂φ) (k × ∇φ) · ê_ω

This is a key prediction of HQIV: vorticity is amplified rather than
decaying as in standard cosmology (where ω ∝ a⁻²).

The source term comes from horizon gradients ∇φ, which are sourced by
density perturbations in the scalar sector.

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


# Load default parameters
_config = load_config()
_chi = float(_config.get('hqiv_parameters', 'chi', fallback='0.172'))
_f_min = float(_config.get('hqiv_parameters', 'f_min', fallback='0.01'))


def compute_df_dphi(alpha, phi, chi=_chi, f_min=_f_min):
    """
    Compute ∂f/∂φ for the vorticity source term.
    
    Paper Reference: paper/main.tex, Section 6
        The vorticity source is proportional to ∂f/∂φ
        
    For the thermodynamic form:
        f = α / (α + χ c φ / 6)
        ∂f/∂φ = -α χ c / 6 / (α + χ c φ / 6)²
               = -χ c / 6 × f × (1 - f) / α
    
    Parameters
    ----------
    alpha : float or array
        Local acceleration magnitude [m/s²]
    phi : float or array
        Horizon field [m/s²]
    chi : float
        Scaling factor
    f_min : float
        Thermodynamic floor
        
    Returns
    -------
    df_dphi : float or array
        Derivative of inertia factor with respect to φ
        
    Notes
    -----
    This derivative is negative: increasing φ (smaller horizon) reduces f.
    The magnitude determines the strength of vorticity generation.
    """
    alpha = np.asarray(alpha)
    phi = np.asarray(phi)
    
    # Compute f first
    a_min = chi * c * phi
    denom = alpha + a_min / 6.0
    denom = np.maximum(denom, 1e-100)
    f = alpha / denom
    f = np.maximum(f, f_min)
    
    # Compute ∂f/∂φ
    # ∂f/∂φ = -α χ c / 6 / (α + χ c φ / 6)²
    #       = -χ c / 6 × f × (1 - f) / α  (alternative form)
    
    df_dphi = -chi * c / 6.0 * f * (1.0 - f) / (alpha + 1e-100)
    
    return df_dphi


def vorticity_source_term(velocity_field, phi, grad_phi, alpha, 
                           chi=_chi, f_min=_f_min, dx=1.0):
    """
    Compute the vorticity source term on the grid.
    
    Paper Reference: paper/main.tex, Section 6
        ∂ω/∂t + (v·∇)ω = (∂f/∂φ) (k × ∇φ) · ê_ω
        
    In the implementation, we compute:
        S_ω = (∂f/∂φ) (k × ∇φ)
        
    where k is related to the local vorticity direction.
    
    Parameters
    ----------
    velocity_field : array
        Velocity field [3, N, N, N]
    phi : array
        Horizon field [N, N, N]
    grad_phi : array
        Gradient of φ [3, N, N, N]
    alpha : array
        Local acceleration magnitude [N, N, N]
    chi : float
        Scaling factor
    f_min : float
        Thermodynamic floor
    dx : float
        Grid spacing [Mpc]
        
    Returns
    -------
    source : array
        Vorticity source term [3, N, N, N]
        
    Notes
    -----
    The source generates vorticity at locations where:
    1. There are gradients in the horizon field (∇φ ≠ 0)
    2. The inertia factor is transitioning (f between f_min and 1)
    
    This happens preferentially at BAO scales during recombination
    and at filament scales during structure formation.
    """
    N = phi.shape[0]
    
    # Compute ∂f/∂φ
    df_dphi = compute_df_dphi(alpha, phi, chi, f_min)
    
    # Compute current vorticity to get direction
    omega = compute_vorticity_field(velocity_field, dx)
    
    # Vorticity magnitude and direction
    omega_mag = np.sqrt(np.sum(omega**2, axis=0)) + 1e-100
    e_omega = omega / omega_mag
    
    # Compute k direction (perpendicular to vorticity)
    # In Fourier space, k is the wavevector
    # For the source, we use the local gradient of velocity
    k_dir = compute_local_k_direction(velocity_field, dx)
    
    # Source: (∂f/∂φ) (k × ∇φ) · ê_ω
    # This is a scalar that multiplies the vorticity direction
    
    # k × ∇φ
    k_cross_grad_phi = np.cross(k_dir, grad_phi, axis=0)
    
    # Dot with vorticity direction
    source_amplitude = df_dphi * np.sum(k_cross_grad_phi * e_omega, axis=0)
    
    # Full source term (in direction of current vorticity)
    source = source_amplitude * e_omega
    
    return source


def compute_local_k_direction(velocity_field, dx):
    """
    Compute the local 'k' direction from the velocity field.
    
    This approximates the local wavevector direction for the vorticity source.
    
    Parameters
    ----------
    velocity_field : array
        Velocity field [3, N, N, N]
    dx : float
        Grid spacing
        
    Returns
    -------
    k_dir : array
        Local k direction [3, N, N, N] (unit vector)
    """
    N = velocity_field.shape[1]
    
    # Compute velocity gradient tensor
    # ∂v_i/∂x_j
    
    # k-vectors for Fourier differentiation
    kx = fftfreq(N, d=dx) * 2 * np.pi
    ky = fftfreq(N, d=dx) * 2 * np.pi
    kz = fftfreq(N, d=dx) * 2 * np.pi
    
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    
    # Transform velocities to Fourier space
    vx_k = fftn(velocity_field[0])
    vy_k = fftn(velocity_field[1])
    vz_k = fftn(velocity_field[2])
    
    # Velocity gradients
    dvx_dx = np.real(ifftn(1j * KX * vx_k))
    dvx_dy = np.real(ifftn(1j * KY * vx_k))
    dvx_dz = np.real(ifftn(1j * KZ * vx_k))
    
    dvy_dx = np.real(ifftn(1j * KX * vy_k))
    dvy_dy = np.real(ifftn(1j * KY * vy_k))
    dvy_dz = np.real(ifftn(1j * KZ * vy_k))
    
    dvz_dx = np.real(ifftn(1j * KX * vz_k))
    dvz_dy = np.real(ifftn(1j * KY * vz_k))
    dvz_dz = np.real(ifftn(1j * KZ * vz_k))
    
    # Strain rate tensor (symmetric part)
    # S_ij = (∂v_i/∂x_j + ∂v_j/∂x_i) / 2
    
    # Principal axes of strain give the local 'k' direction
    # For simplicity, use the direction of maximum velocity gradient
    
    grad_v_mag = np.sqrt(
        dvx_dx**2 + dvx_dy**2 + dvx_dz**2 +
        dvy_dx**2 + dvy_dy**2 + dvy_dz**2 +
        dvz_dx**2 + dvz_dy**2 + dvz_dz**2
    )
    
    # Use velocity direction as proxy for k direction
    v_mag = np.sqrt(np.sum(velocity_field**2, axis=0)) + 1e-100
    k_dir = velocity_field / v_mag
    
    return k_dir


def compute_vorticity_field(velocity_field, dx):
    """
    Compute vorticity field ω = ∇ × v.
    
    Parameters
    ----------
    velocity_field : array
        Velocity field [3, N, N, N]
    dx : float
        Grid spacing
        
    Returns
    -------
    omega : array
        Vorticity field [3, N, N, N]
        
    Notes
    -----
    In Fourier space:
        ω_k = i k × v_k
    """
    N = velocity_field.shape[1]
    
    # k-vectors
    kx = fftfreq(N, d=dx) * 2 * np.pi
    ky = fftfreq(N, d=dx) * 2 * np.pi
    kz = fftfreq(N, d=dx) * 2 * np.pi
    
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    
    # Transform velocities to Fourier space
    vx_k = fftn(velocity_field[0])
    vy_k = fftn(velocity_field[1])
    vz_k = fftn(velocity_field[2])
    
    # Vorticity in k-space: ω_k = i k × v_k
    omega_x_k = 1j * (KY * vz_k - KZ * vy_k)
    omega_y_k = 1j * (KZ * vx_k - KX * vz_k)
    omega_z_k = 1j * (KX * vy_k - KY * vx_k)
    
    # Transform back to real space
    omega_x = np.real(ifftn(omega_x_k))
    omega_y = np.real(ifftn(omega_y_k))
    omega_z = np.real(ifftn(omega_z_k))
    
    return np.stack([omega_x, omega_y, omega_z])


def compute_vorticity_evolution(omega, source, H_a, a, dt):
    """
    Evolve vorticity field by one timestep.
    
    Paper equation:
        ∂ω/∂t + (v·∇)ω = source
        
    In scale factor time:
        dω/da = source / (a H) - 2 ω / a
        
    The -2 ω / a term is the standard decay in expanding universe.
    The source term can counteract this decay.
    
    Parameters
    ----------
    omega : array
        Current vorticity field [3, N, N, N]
    source : array
        Vorticity source term [3, N, N, N]
    H_a : float
        Hubble parameter [1/s]
    a : float
        Scale factor
    dt : float
        Time step [s]
        
    Returns
    -------
    omega_new : array
        Updated vorticity field
    """
    # Standard decay: dω/dt = -2 H ω (in expanding universe)
    decay_rate = 2.0 * H_a
    
    # Source contribution
    # The source can be positive (amplification) or negative (damping)
    
    # Euler step
    # dω/dt = source - 2 H ω
    domega_dt = source - decay_rate * omega
    
    omega_new = omega + domega_dt * dt
    
    return omega_new


def vorticity_growth_exponent(omega_history, a_history):
    """
    Compute the effective vorticity growth exponent.
    
    Paper prediction: ω ∝ a^β for some effective exponent β.
    In standard cosmology: β = -2 (decay).
    In HQIV: β can be positive (amplification).
    
    Parameters
    ----------
    omega_history : list
        List of vorticity magnitude arrays at different times
    a_history : array
        Scale factors corresponding to omega_history
        
    Returns
    -------
    beta_eff : float
        Effective growth exponent
    """
    # Compute mean vorticity magnitude at each time
    omega_mag_mean = np.array([
        np.mean(np.sqrt(np.sum(omega**2, axis=0)))
        for omega in omega_history
    ])
    
    # Fit power law: ω ∝ a^β
    # log(ω) = β log(a) + const
    log_a = np.log(a_history)
    log_omega = np.log(omega_mag_mean + 1e-100)
    
    # Linear fit
    coeffs = np.polyfit(log_a, log_omega, 1)
    beta_eff = coeffs[0]
    
    return beta_eff


class VorticitySolver:
    """
    Class to manage vorticity computation throughout the simulation.
    
    Parameters
    ----------
    H0 : float
        Hubble constant [km/s/Mpc]
    chi : float
        Scaling factor
    f_min : float
        Thermodynamic floor
    N_grid : int
        Grid size
    box_size : float
        Box size [Mpc]
    """
    
    def __init__(self, H0=73.2, chi=_chi, f_min=_f_min,
                 N_grid=256, box_size=100.0):
        self.H0 = H0 * 1000.0 / Mpc_m  # Convert to 1/s
        self.chi = chi
        self.f_min = f_min
        self.N_grid = N_grid
        self.box_size = box_size
        self.dx = box_size / N_grid
        
        # Storage
        self.omega = None
        self.source = None
        self.omega_history = []
        self.a_history = []
        
    def compute_vorticity(self, velocity_field):
        """Compute vorticity from velocity field."""
        self.omega = compute_vorticity_field(velocity_field, self.dx)
        return self.omega
    
    def compute_source(self, velocity_field, phi, grad_phi, alpha):
        """Compute vorticity source term."""
        self.source = vorticity_source_term(
            velocity_field, phi, grad_phi, alpha,
            self.chi, self.f_min, self.dx
        )
        return self.source
    
    def evolve(self, H_a, a, dt):
        """Evolve vorticity by one timestep."""
        if self.omega is None or self.source is None:
            raise RuntimeError("Call compute_vorticity and compute_source first")
        
        self.omega = compute_vorticity_evolution(
            self.omega, self.source, H_a, a, dt
        )
        
        # Store history
        self.omega_history.append(self.omega.copy())
        self.a_history.append(a)
        
        return self.omega
    
    def get_growth_exponent(self):
        """Compute effective growth exponent from history."""
        if len(self.omega_history) < 2:
            return None
        
        return vorticity_growth_exponent(
            self.omega_history, 
            np.array(self.a_history)
        )
    
    def apply_to_velocities(self, velocities, positions, dt):
        """
        Apply vorticity change to particle velocities.
        
        The vorticity change induces a velocity change:
            Δv ≈ Δω × r
        
        Parameters
        ----------
        velocities : array
            Particle velocities [N_part, 3]
        positions : array
            Particle positions [N_part, 3]
        dt : float
            Time step
            
        Returns
        -------
        velocities_new : array
            Updated velocities
        """
        if self.omega is None:
            return velocities
        
        # Interpolate vorticity to particle positions
        omega_p = self._interpolate_to_particles(self.omega, positions)
        
        # Position relative to box center
        center = self.box_size / 2
        r_rel = positions - center
        
        # Velocity change from vorticity: Δv = ω × r × dt
        # This adds rotational motion
        dv = np.cross(omega_p, r_rel) * dt * 0.1  # Scale factor for stability
        
        velocities_new = velocities + dv
        
        return velocities_new
    
    def _interpolate_to_particles(self, field, positions):
        """Interpolate field to particle positions using CIC."""
        values = np.zeros((len(positions), 3))
        
        for i, pos in enumerate(positions):
            ix = int(np.floor(pos[0] / self.dx)) % self.N_grid
            iy = int(np.floor(pos[1] / self.dx)) % self.N_grid
            iz = int(np.floor(pos[2] / self.dx)) % self.N_grid
            
            ux = pos[0] / self.dx - np.floor(pos[0] / self.dx)
            uy = pos[1] / self.dx - np.floor(pos[1] / self.dx)
            uz = pos[2] / self.dx - np.floor(pos[2] / self.dx)
            
            for ddx in [0, 1]:
                for ddy in [0, 1]:
                    for ddz in [0, 1]:
                        w = ((1-ux) if ddx==0 else ux) * \
                            ((1-uy) if ddy==0 else uy) * \
                            ((1-uz) if ddz==0 else uz)
                        for d in range(3):
                            values[i, d] += w * field[d, (ix+ddx)%self.N_grid, (iy+ddy)%self.N_grid, (iz+ddz)%self.N_grid]
        
        return values


# =============================================================================
# Testing
# =============================================================================

def test_vorticity():
    """Test vorticity computation."""
    import matplotlib.pyplot as plt
    
    print("Testing vorticity computation...")
    
    # Create test velocity field with rotation
    N = 64
    box_size = 100.0
    dx = box_size / N
    
    x = np.linspace(0, box_size, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    center = box_size / 2
    
    # Rotational velocity field
    velocity = np.zeros((3, N, N, N))
    omega_test = 0.1  # Test vorticity amplitude
    
    # v = ω × r (solid body rotation around z-axis)
    velocity[0] = -omega_test * (Y - center)
    velocity[1] = omega_test * (X - center)
    velocity[2] = 0
    
    # Compute vorticity
    omega = compute_vorticity_field(velocity, dx)
    omega_mag = np.sqrt(np.sum(omega**2, axis=0))
    
    print(f"Expected ω_z = {2*omega_test:.4f}")
    print(f"Computed ω_z (mean) = {np.mean(omega[2]):.4f}")
    print(f"ω_z (std) = {np.std(omega[2]):.4f}")
    
    # Test source term
    phi = np.ones((N, N, N)) * 1e-10  # Constant φ
    grad_phi = np.zeros((3, N, N, N))
    grad_phi[0] = 1e-12  # Small gradient
    
    alpha = np.ones((N, N, N)) * 1e-10
    
    source = vorticity_source_term(velocity, phi, grad_phi, alpha, dx=dx)
    source_mag = np.sqrt(np.sum(source**2, axis=0))
    
    print(f"Source magnitude (mean) = {np.mean(source_mag):.2e}")
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    ax = axes[0]
    im = ax.imshow(velocity[0, :, :, N//2], origin='lower', extent=[0, box_size, 0, box_size])
    ax.set_title('v_x (slice)')
    plt.colorbar(im, ax=ax)
    
    ax = axes[1]
    im = ax.imshow(omega[2, :, :, N//2], origin='lower', extent=[0, box_size, 0, box_size])
    ax.set_title('ω_z (slice)')
    plt.colorbar(im, ax=ax)
    
    ax = axes[2]
    im = ax.imshow(source_mag[:, :, N//2], origin='lower', extent=[0, box_size, 0, box_size])
    ax.set_title('Source magnitude (slice)')
    plt.colorbar(im, ax=ax)
    
    plt.tight_layout()
    plt.savefig('vorticity_test.png', dpi=150)
    plt.close()
    print("Saved vorticity_test.png")
    
    return omega


if __name__ == "__main__":
    test_vorticity()