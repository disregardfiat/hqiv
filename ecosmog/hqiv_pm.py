"""
HQIV Particle-Mesh N-body Simulator

A Python implementation of the HQIV cosmology for structure formation tests.
Uses PM (Particle-Mesh) method with:
- Modified gravitational coupling G_eff(a,k)
- Horizon term in Poisson equation
- Scale-dependent inertia reduction

Author: HQIV Team
"""

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq, rfftn, irfftn
import os
import sys

# Add sandbox to path for background solver
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'sandbox'))
from hqiv_background import run_background, H0, Omega_m_vis, Omega_r0, beta, c, Mpc, G0

try:
    import matplotlib.pyplot as plt
    _has_plt = True
except ImportError:
    _has_plt = False


class HQIVSimulator:
    """
    Particle-Mesh N-body simulator for HQIV cosmology.
    
    Uses comoving coordinates where:
    - x_comoving = x_physical / a
    - v = a * dx_comoving/dt = peculiar velocity
    - ∇²Φ = 4πG ρ_mean a δ (Poisson equation in comoving coords)
    
    Equations of motion:
    - dx/dt = v / a
    - dv/dt = -H v - ∇Φ / a
    
    Vorticity equation (from paper):
    - ∂ω/∂t + (v·∇)ω = β H (ω · ê_Θ)
    - This amplifies vorticity aligned with horizon direction
    """
    
    def __init__(self, N_particles=64**3, L_box=100.0, N_grid=64, 
                 a_start=0.02, a_end=1.0, n_steps=100):
        """
        Initialize the simulator.
        """
        self.N_particles = N_particles
        self.L_box = L_box
        self.N_grid = N_grid
        self.a_start = a_start
        self.a_end = a_end
        self.n_steps = n_steps
        
        # Grid spacing in comoving coordinates
        self.dx = L_box / N_grid
        
        # Particle mass (total mass = 1 in code units)
        self.particle_mass = 1.0 / N_particles
        
        # Load background H(a) table
        print("Loading background H(a) table...")
        self.a_bg, self.H_bg = run_background(a_min=1e-6, a_max=1.0, n_pts=2000)
        
        # Initialize particles
        self._init_particles()
        
        # Precompute k-space grid
        self._init_kspace()
        
        print(f"HQIV PM Simulator initialized:")
        print(f"  N_particles = {N_particles}")
        print(f"  L_box = {L_box} Mpc/h")
        print(f"  N_grid = {N_grid}")
        print(f"  a_start = {a_start}, a_end = {a_end}")
    
    def _init_particles(self):
        """Initialize particles with Zel'dovich approximation."""
        n_per_dim = int(np.round(self.N_particles ** (1/3)))
        
        # Create uniform grid in comoving coordinates
        x = np.linspace(0, self.L_box, n_per_dim, endpoint=False)
        y = np.linspace(0, self.L_box, n_per_dim, endpoint=False)
        z = np.linspace(0, self.L_box, n_per_dim, endpoint=False)
        
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        self.pos = np.zeros((self.N_particles, 3))
        self.pos[:, 0] = X.flatten()[:self.N_particles]
        self.pos[:, 1] = Y.flatten()[:self.N_particles]
        self.pos[:, 2] = Z.flatten()[:self.N_particles]
        
        # Zel'dovich approximation for initial conditions
        np.random.seed(42)
        
        # Generate Gaussian random field with CDM-like power spectrum
        kx = fftfreq(n_per_dim, d=1.0/n_per_dim) * 2 * np.pi / self.L_box
        ky = fftfreq(n_per_dim, d=1.0/n_per_dim) * 2 * np.pi / self.L_box
        kz = fftfreq(n_per_dim, d=1.0/n_per_dim) * 2 * np.pi / self.L_box
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        K = np.sqrt(KX**2 + KY**2 + KZ**2)
        K[0, 0, 0] = 1.0
        
        # CDM-like power spectrum: P(k) ~ k^n with transfer function
        # Use a reasonable amplitude for visible structure formation
        n_index = -1.0  # roughly scale-invariant
        k_nyq = np.pi * n_per_dim / self.L_box
        
        # Power spectrum with small-scale cutoff
        # Amplitude scaled to give ~10% density fluctuations at start
        P_k = 0.1 * K**(n_index) * np.exp(-(K/k_nyq)**2)
        P_k[0, 0, 0] = 0
        
        # Generate complex Gaussian random field
        phases = np.random.uniform(0, 2*np.pi, (n_per_dim, n_per_dim, n_per_dim))
        amplitudes = np.random.rayleigh(scale=1.0, size=(n_per_dim, n_per_dim, n_per_dim))
        
        # Complex field in k-space
        delta_k = amplitudes * np.sqrt(P_k) * np.exp(1j * phases)
        
        # Displacement field: ψ_k = -i k δ_k / k² (from Zel'dovich)
        psi_x_k = -1j * KX * delta_k / (K**2)
        psi_y_k = -1j * KY * delta_k / (K**2)
        psi_z_k = -1j * KZ * delta_k / (K**2)
        
        # Zero out k=0 mode
        psi_x_k[0, 0, 0] = 0
        psi_y_k[0, 0, 0] = 0
        psi_z_k[0, 0, 0] = 0
        
        # Transform to real space
        psi_x = np.real(ifftn(psi_x_k))
        psi_y = np.real(ifftn(psi_y_k))
        psi_z = np.real(ifftn(psi_z_k))
        
        # Apply Zel'dovich displacement
        # x = x_grid + D(a) * ψ
        # where D(a) is the linear growth factor
        D_i = self.a_start  # In matter domination, D ~ a
        
        self.pos[:, 0] += D_i * psi_x.flatten()[:self.N_particles]
        self.pos[:, 1] += D_i * psi_y.flatten()[:self.N_particles]
        self.pos[:, 2] += D_i * psi_z.flatten()[:self.N_particles]
        
        # Apply periodic boundary conditions
        self.pos = self.pos % self.L_box
        
        # Initial peculiar velocity: v = a * H * dD/da * ψ = a * H * ψ (for D ~ a)
        H_a = self.H_of_a(self.a_start)
        f_growth = 1.0  # d ln D / d ln a = 1 for matter domination
        
        self.vel = np.zeros((self.N_particles, 3))
        self.vel[:, 0] = self.a_start * H_a * f_growth * psi_x.flatten()[:self.N_particles]
        self.vel[:, 1] = self.a_start * H_a * f_growth * psi_y.flatten()[:self.N_particles]
        self.vel[:, 2] = self.a_start * H_a * f_growth * psi_z.flatten()[:self.N_particles]
        
        print(f"  Initial displacement RMS: {np.std(psi_x):.2e}")
        print(f"  Initial velocity RMS: {np.std(self.vel):.2e}")
    
    def _init_kspace(self):
        """Initialize k-space grid for FFT-based Poisson solver."""
        # k values in comoving units (h/Mpc)
        self.kx = fftfreq(self.N_grid, d=self.dx) * 2 * np.pi
        self.ky = fftfreq(self.N_grid, d=self.dx) * 2 * np.pi
        self.kz = fftfreq(self.N_grid, d=self.dx) * 2 * np.pi
        
        KX, KY, KZ = np.meshgrid(self.kx, self.ky, self.kz, indexing='ij')
        
        # k^2 in (h/Mpc)^2
        self.k2 = KX**2 + KY**2 + KZ**2
        self.k2[0, 0, 0] = 1.0  # Avoid division by zero
        
        # k magnitude
        self.k_mag = np.sqrt(self.k2)
    
    def H_of_a(self, a):
        """Interpolate H(a) from background table."""
        return np.interp(a, self.a_bg, self.H_bg)
    
    def G_eff(self, a, k):
        """
        Effective gravitational coupling in HQIV.
        
        G_eff/G0 = (H(a)/H0)^α * horizon_factor(k, a)
        """
        H_a = self.H_of_a(a)
        
        # G ratio from varying G (α = 0.6 from QI)
        alpha = 0.6
        G_ratio = (H_a / H0) ** alpha
        
        # Horizon cutoff
        Theta = 2 * c / H_a  # horizon in meters
        Theta_Mpc = Theta / Mpc  # in Mpc
        k_horizon = 2 * np.pi / Theta_Mpc  # horizon wavenumber
        
        # Smooth suppression for k < k_horizon
        horizon_factor = np.ones_like(k)
        mask = k < k_horizon
        if np.any(mask):
            horizon_factor[mask] = (k[mask] / k_horizon)**0.5
        
        return G_ratio * horizon_factor
    
    def compute_density(self):
        """Compute density contrast δ on grid using CIC."""
        rho = np.zeros((self.N_grid, self.N_grid, self.N_grid))
        
        # CIC assignment
        for i in range(self.N_particles):
            x, y, z = self.pos[i]
            
            # Cell indices and weights
            ix = int(np.floor(x / self.dx)) % self.N_grid
            iy = int(np.floor(y / self.dx)) % self.N_grid
            iz = int(np.floor(z / self.dx)) % self.N_grid
            
            ux = x / self.dx - np.floor(x / self.dx)
            uy = y / self.dx - np.floor(y / self.dx)
            uz = z / self.dx - np.floor(z / self.dx)
            
            # Distribute to 8 cells
            for dx in [0, 1]:
                for dy in [0, 1]:
                    for dz in [0, 1]:
                        w = ((1-ux) if dx==0 else ux) * \
                            ((1-uy) if dy==0 else uy) * \
                            ((1-uz) if dz==0 else uz)
                        rho[(ix+dx)%self.N_grid, (iy+dy)%self.N_grid, (iz+dz)%self.N_grid] += w
        
        # Normalize to get density contrast
        rho_mean = self.N_particles / self.N_grid**3
        delta = rho / rho_mean - 1.0
        
        return delta
    
    def solve_poisson(self, delta, a):
        """
        Solve Poisson equation for gravitational potential.
        
        In comoving coordinates with our normalization:
        ∇²Φ = 4πG ρ_mean a³ δ
        
        The force is then: F = -∇Φ
        And acceleration: a = F / a (peculiar acceleration)
        
        In k-space:
        Φ_k = -4πG ρ_mean a³ δ_k / k²
        
        With our code units (total mass = 1, L_box = 1 in code units):
        4πG ρ_mean = 4πG * (1/L³) = 3 Ω_m H0² / (2 a³)
        
        So: Φ_k = -(3/2) Ω_m H0² a³ δ_k / (a³ k²) = -(3/2) Ω_m δ_k / k²
        
        The force is: F_k = -i k Φ_k
        """
        delta_k = fftn(delta)
        
        # Scale-dependent G_eff
        G_eff_k = self.G_eff(a, self.k_mag)
        
        # Gravitational potential in k-space
        # Φ_k = -(3/2) Ω_m G_eff δ_k / k²
        phi_k = np.zeros_like(delta_k, dtype=complex)
        nonzero = self.k2 > 0
        
        # Simple, stable normalization
        phi_k[nonzero] = -1.5 * Omega_m_vis * G_eff_k[nonzero] * delta_k[nonzero] / self.k2[nonzero]
        phi_k[0, 0, 0] = 0
        
        return phi_k
    
    def compute_acceleration(self, phi_k, a):
        """
        Compute gravitational acceleration in comoving coordinates.
        
        g = -∇Φ / a  (peculiar acceleration)
        
        In k-space: g_k = -i k Φ_k / a
        """
        KX, KY, KZ = np.meshgrid(self.kx, self.ky, self.kz, indexing='ij')
        
        # Acceleration in k-space
        gx_k = -1j * KX * phi_k / a
        gy_k = -1j * KY * phi_k / a
        gz_k = -1j * KZ * phi_k / a
        
        # Transform to real space
        gx = np.real(ifftn(gx_k))
        gy = np.real(ifftn(gy_k))
        gz = np.real(ifftn(gz_k))
        
        return gx, gy, gz
    
    def interpolate_acceleration(self, gx, gy, gz, pos):
        """Interpolate acceleration field to particle positions using CIC."""
        acc = np.zeros((len(pos), 3))
        
        for i in range(len(pos)):
            x, y, z = pos[i]
            
            ix = int(np.floor(x / self.dx)) % self.N_grid
            iy = int(np.floor(y / self.dx)) % self.N_grid
            iz = int(np.floor(z / self.dx)) % self.N_grid
            
            ux = x / self.dx - np.floor(x / self.dx)
            uy = y / self.dx - np.floor(y / self.dx)
            uz = z / self.dx - np.floor(z / self.dx)
            
            for dx in [0, 1]:
                for dy in [0, 1]:
                    for dz in [0, 1]:
                        w = ((1-ux) if dx==0 else ux) * \
                            ((1-uy) if dy==0 else uy) * \
                            ((1-uz) if dz==0 else uz)
                        acc[i, 0] += w * gx[(ix+dx)%self.N_grid, (iy+dy)%self.N_grid, (iz+dz)%self.N_grid]
                        acc[i, 1] += w * gy[(ix+dx)%self.N_grid, (iy+dy)%self.N_grid, (iz+dz)%self.N_grid]
                        acc[i, 2] += w * gz[(ix+dx)%self.N_grid, (iy+dy)%self.N_grid, (iz+dz)%self.N_grid]
        
        return acc
    
    def inertia_factor(self, a, acc_mag):
        """
        Compute inertia reduction factor from QI.
        
        m_i/m_g = 1 - a_min/|a_local|
        
        where a_min = β c H(a)
        """
        H_a = self.H_of_a(a)
        a_min = beta * c * H_a
        
        # Inertia factor
        f = 1.0 - a_min / (acc_mag + 1e-100)
        
        # Floor to prevent instability
        return np.maximum(f, 0.35)
    
    def step(self, a, da):
        """
        Advance simulation by one time step using kick-drift-kick leapfrog.
        
        We work in dimensionless code units where:
        - Length is in units of L_box
        - Time is in units of 1/H0
        - Velocity is in units of L_box * H0
        
        Equations of motion in scale factor time:
        - dx/da = v / (H' a²)  where H' = H/H0
        - dv/da = -v/(H' a) + g'/(H' a)
        
        where g' is the dimensionless acceleration.
        """
        H_a = self.H_of_a(a)
        H_prime = H_a / H0  # dimensionless Hubble parameter
        
        # Compute density and potential
        delta = self.compute_density()
        phi_k = self.solve_poisson(delta, a)
        
        # Compute acceleration field
        gx, gy, gz = self.compute_acceleration(phi_k, a)
        
        # Interpolate to particles
        acc = self.interpolate_acceleration(gx, gy, gz, self.pos)
        
        # Apply inertia modification
        acc_mag = np.sqrt(np.sum(acc**2, axis=1))
        f_inertia = self.inertia_factor(a, acc_mag)
        acc *= (1.0 / f_inertia)[:, np.newaxis]
        
        # Limit acceleration to prevent instability
        acc_max = 1e3
        acc_mag = np.sqrt(np.sum(acc**2, axis=1))
        too_large = acc_mag > acc_max
        if np.any(too_large):
            acc[too_large] *= acc_max / acc_mag[too_large, np.newaxis]
        
        # Kick-drift-kick leapfrog in dimensionless units
        # First half kick
        self.vel += 0.5 * da * (-self.vel / a + acc) / H_prime
        
        # Limit velocity
        vel_max = 1e2
        vel_mag = np.sqrt(np.sum(self.vel**2, axis=1))
        too_fast = vel_mag > vel_max
        if np.any(too_fast):
            self.vel[too_fast] *= vel_max / vel_mag[too_fast, np.newaxis]
        
        # Drift
        self.pos += da * self.vel / (H_prime * a * a)
        self.pos = self.pos % self.L_box
        
        # Recompute acceleration at new position
        a_new = a + da
        H_new = self.H_of_a(a_new)
        H_prime_new = H_new / H0
        
        delta = self.compute_density()
        phi_k = self.solve_poisson(delta, a_new)
        gx, gy, gz = self.compute_acceleration(phi_k, a_new)
        acc = self.interpolate_acceleration(gx, gy, gz, self.pos)
        
        acc_mag = np.sqrt(np.sum(acc**2, axis=1))
        f_inertia = self.inertia_factor(a_new, acc_mag)
        acc *= (1.0 / f_inertia)[:, np.newaxis]
        
        # Limit acceleration
        acc_mag = np.sqrt(np.sum(acc**2, axis=1))
        too_large = acc_mag > acc_max
        if np.any(too_large):
            acc[too_large] *= acc_max / acc_mag[too_large, np.newaxis]
        
        # Second half kick
        self.vel += 0.5 * da * (-self.vel / a_new + acc) / H_prime_new
        
        # Limit velocity again
        vel_mag = np.sqrt(np.sum(self.vel**2, axis=1))
        too_fast = vel_mag > vel_max
        if np.any(too_fast):
            self.vel[too_fast] *= vel_max / vel_mag[too_fast, np.newaxis]
    
    def compute_velocity_field(self):
        """Compute velocity field on grid using CIC from particle velocities."""
        vx = np.zeros((self.N_grid, self.N_grid, self.N_grid))
        vy = np.zeros((self.N_grid, self.N_grid, self.N_grid))
        vz = np.zeros((self.N_grid, self.N_grid, self.N_grid))
        weights = np.zeros((self.N_grid, self.N_grid, self.N_grid))
        
        # CIC assignment for velocities
        for i in range(self.N_particles):
            x, y, z = self.pos[i]
            vx_p, vy_p, vz_p = self.vel[i]
            
            ix = int(np.floor(x / self.dx)) % self.N_grid
            iy = int(np.floor(y / self.dx)) % self.N_grid
            iz = int(np.floor(z / self.dx)) % self.N_grid
            
            ux = x / self.dx - np.floor(x / self.dx)
            uy = y / self.dx - np.floor(y / self.dx)
            uz = z / self.dx - np.floor(z / self.dx)
            
            for dx in [0, 1]:
                for dy in [0, 1]:
                    for dz in [0, 1]:
                        w = ((1-ux) if dx==0 else ux) * \
                            ((1-uy) if dy==0 else uy) * \
                            ((1-uz) if dz==0 else uz)
                        ix2 = (ix+dx) % self.N_grid
                        iy2 = (iy+dy) % self.N_grid
                        iz2 = (iz+dz) % self.N_grid
                        vx[ix2, iy2, iz2] += w * vx_p
                        vy[ix2, iy2, iz2] += w * vy_p
                        vz[ix2, iy2, iz2] += w * vz_p
                        weights[ix2, iy2, iz2] += w
        
        # Normalize by weights
        mask = weights > 0
        vx[mask] /= weights[mask]
        vy[mask] /= weights[mask]
        vz[mask] /= weights[mask]
        
        return vx, vy, vz
    
    def compute_vorticity_field(self, vx, vy, vz):
        """
        Compute vorticity field ω = ∇ × v.
        
        In Fourier space:
        ω_k = i k × v_k
        
        Components:
        ω_x = ∂v_z/∂y - ∂v_y/∂z
        ω_y = ∂v_x/∂z - ∂v_z/∂x
        ω_z = ∂v_y/∂x - ∂v_x/∂y
        """
        # Transform velocities to k-space
        vx_k = fftn(vx)
        vy_k = fftn(vy)
        vz_k = fftn(vz)
        
        # k vectors
        KX, KY, KZ = np.meshgrid(self.kx, self.ky, self.kz, indexing='ij')
        
        # Vorticity in k-space: ω_k = i k × v_k
        # ω_x = i (k_y * v_z_k - k_z * v_y_k)
        # ω_y = i (k_z * v_x_k - k_x * v_z_k)
        # ω_z = i (k_x * v_y_k - k_y * v_x_k)
        omega_x_k = 1j * (KY * vz_k - KZ * vy_k)
        omega_y_k = 1j * (KZ * vx_k - KX * vz_k)
        omega_z_k = 1j * (KX * vy_k - KY * vx_k)
        
        # Transform back to real space
        omega_x = np.real(ifftn(omega_x_k))
        omega_y = np.real(ifftn(omega_y_k))
        omega_z = np.real(ifftn(omega_z_k))
        
        return omega_x, omega_y, omega_z
    
    def compute_horizon_direction(self, pos):
        """
        Compute horizon direction unit vector ê_Θ for each particle.
        
        The horizon direction points radially outward from the observer
        toward the particle's position. In a homogeneous universe,
        this is the direction to the particle's causal horizon.
        
        For simplicity, we use the position direction normalized.
        """
        # Normalize position to get horizon direction
        r = np.sqrt(np.sum(pos**2, axis=1, keepdims=True)) + 1e-100
        e_theta = pos / r
        
        return e_theta
    
    def apply_vorticity_amplification(self, a, da):
        """
        Apply vorticity amplification from horizon coupling.
        
        From paper: ∂ω/∂t + (v·∇)ω = β H (ω · ê_Θ)
        
        In scale factor time:
        ∂ω/∂a = β (ω · ê_Θ) / (a H)
        
        This amplifies vorticity components aligned with the horizon direction.
        """
        H_a = self.H_of_a(a)
        
        # Compute velocity field
        vx, vy, vz = self.compute_velocity_field()
        
        # Compute vorticity field
        omega_x, omega_y, omega_z = self.compute_vorticity_field(vx, vy, vz)
        
        # Interpolate vorticity to particle positions
        omega_p = np.zeros((self.N_particles, 3))
        for i in range(self.N_particles):
            x, y, z = self.pos[i]
            
            ix = int(np.floor(x / self.dx)) % self.N_grid
            iy = int(np.floor(y / self.dx)) % self.N_grid
            iz = int(np.floor(z / self.dx)) % self.N_grid
            
            ux = x / self.dx - np.floor(x / self.dx)
            uy = y / self.dx - np.floor(y / self.dx)
            uz = z / self.dx - np.floor(z / self.dx)
            
            for dx in [0, 1]:
                for dy in [0, 1]:
                    for dz in [0, 1]:
                        w = ((1-ux) if dx==0 else ux) * \
                            ((1-uy) if dy==0 else uy) * \
                            ((1-uz) if dz==0 else uz)
                        ix2 = (ix+dx) % self.N_grid
                        iy2 = (iy+dy) % self.N_grid
                        iz2 = (iz+dz) % self.N_grid
                        omega_p[i, 0] += w * omega_x[ix2, iy2, iz2]
                        omega_p[i, 1] += w * omega_y[ix2, iy2, iz2]
                        omega_p[i, 2] += w * omega_z[ix2, iy2, iz2]
        
        # Compute horizon direction for each particle
        e_theta = self.compute_horizon_direction(self.pos)
        
        # Compute ω · ê_Θ (projection of vorticity onto horizon direction)
        omega_dot_e = np.sum(omega_p * e_theta, axis=1)
        
        # Amplification factor: Δω = β * (ω · ê_Θ) * ê_Θ * da / (a * H)
        # This adds vorticity in the horizon direction
        amplification = beta * omega_dot_e[:, np.newaxis] * e_theta * da / (a * H_a)
        
        # Apply amplification to velocities
        # The vorticity change translates to a velocity change
        # For small changes: Δv ≈ Δω × r (where r is position relative to center)
        # Simpler approach: directly modify velocity based on vorticity change
        self.vel += amplification * 0.1  # Scale factor for stability
        
        # Return vorticity statistics
        omega_mag = np.sqrt(np.sum(omega_p**2, axis=1))
        return {
            'omega_rms': np.sqrt(np.mean(omega_mag**2)),
            'omega_max': np.max(omega_mag),
            'omega_dot_e_mean': np.mean(np.abs(omega_dot_e))
        }
    
    def compute_angular_momentum(self):
        """
        Compute angular momentum for each particle.
        
        L = m * r × v (in comoving coordinates)
        
        Returns total angular momentum and per-particle values.
        """
        # Position relative to box center
        center = self.L_box / 2
        r_rel = self.pos - center
        
        # Angular momentum: L = r × v
        L = np.cross(r_rel, self.vel)
        
        # Total angular momentum
        L_total = np.sum(L, axis=0) * self.particle_mass
        
        # Angular momentum magnitude per particle
        L_mag = np.sqrt(np.sum(L**2, axis=1))
        
        return L_total, L_mag
    
    def compute_spin_parameter(self):
        """
        Compute the spin parameter λ for the entire system.
        
        λ = |J| * sqrt(|E|) / (G * M^(5/2))
        
        where J is total angular momentum, E is total energy,
        and M is total mass.
        
        This is a dimensionless measure of rotation.
        """
        # Total angular momentum
        L_total, L_mag = self.compute_angular_momentum()
        J = np.sqrt(np.sum(L_total**2))
        
        # Kinetic energy
        KE = 0.5 * np.sum(np.sum(self.vel**2, axis=1)) * self.particle_mass
        
        # Potential energy (approximate from density field)
        delta = self.compute_density()
        phi_k = self.solve_poisson(delta, 1.0)
        phi = np.real(ifftn(phi_k))
        
        # PE ≈ 0.5 * Σ ρ * Φ (virial theorem)
        # Simplified: use kinetic energy as proxy
        E_total = KE  # Approximation
        
        # Total mass
        M_total = self.N_particles * self.particle_mass
        
        # Spin parameter (simplified)
        if E_total > 0:
            lambda_spin = J * np.sqrt(E_total) / (M_total**(5/2) + 1e-100)
        else:
            lambda_spin = 0
        
        return lambda_spin, J, KE
    
    def compute_power_spectrum(self):
        """Compute matter power spectrum P(k)."""
        delta = self.compute_density()
        delta_k = fftn(delta)
        
        # Power in each k-mode
        power = np.abs(delta_k)**2
        
        # Bin by k magnitude
        k_max = np.max(self.k_mag)
        k_bins = np.linspace(0, k_max, 50)
        k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
        P_k = np.zeros(len(k_centers))
        
        k_flat = self.k_mag.flatten()
        p_flat = power.flatten()
        
        for i in range(len(k_centers)):
            mask = (k_flat >= k_bins[i]) & (k_flat < k_bins[i+1])
            if np.any(mask):
                P_k[i] = np.mean(p_flat[mask])
        
        # Normalize: P(k) = <|δ_k|²> * V
        P_k *= self.L_box**3
        
        return k_centers, P_k
    
    def run(self, save_snapshots=True, snapshot_interval=10, enable_vorticity=True):
        """
        Run the simulation with vorticity amplification.
        
        Parameters:
        -----------
        save_snapshots : bool
            Save power spectrum snapshots
        snapshot_interval : int
            Interval between snapshots
        enable_vorticity : bool
            Enable vorticity amplification from horizon coupling
        """
        print(f"\nRunning HQIV PM simulation...")
        print(f"  From a={self.a_start} to a={self.a_end} in {self.n_steps} steps")
        print(f"  Vorticity amplification: {'ON' if enable_vorticity else 'OFF'}")
        
        a_arr = np.linspace(self.a_start, self.a_end, self.n_steps + 1)
        
        self.power_spectra = []
        self.scale_factors = []
        self.vorticity_stats = []
        self.angular_momentum = []
        self.spin_parameters = []
        
        for i in range(self.n_steps):
            a = a_arr[i]
            da = a_arr[i+1] - a_arr[i]
            
            # Standard PM step
            self.step(a, da)
            
            # Apply vorticity amplification
            if enable_vorticity:
                vort_stats = self.apply_vorticity_amplification(a, da)
            
            if (i + 1) % snapshot_interval == 0:
                # Power spectrum
                k, Pk = self.compute_power_spectrum()
                self.power_spectra.append((k, Pk))
                self.scale_factors.append(a_arr[i+1])
                
                # Vorticity statistics
                if enable_vorticity:
                    self.vorticity_stats.append(vort_stats)
                
                # Angular momentum and spin
                L_total, L_mag = self.compute_angular_momentum()
                lambda_spin, J, KE = self.compute_spin_parameter()
                self.angular_momentum.append(L_total)
                self.spin_parameters.append(lambda_spin)
                
                print(f"  Step {i+1}/{self.n_steps}, a={a_arr[i+1]:.4f}, "
                      f"P(k peak)={np.max(Pk):.2e}, "
                      f"λ={lambda_spin:.2e}, "
                      f"ω_rms={vort_stats['omega_rms']:.2e}" if enable_vorticity else "")
        
        print("Simulation complete!")
        return self.power_spectra, self.scale_factors
    
    def plot_power_spectrum(self, output_file='hqiv_power_spectrum.png'):
        """Plot the power spectrum evolution."""
        if not _has_plt:
            print("matplotlib not available, skipping plot")
            return
        
        if len(self.power_spectra) == 0:
            print("No power spectra to plot. Run simulation first.")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Power spectra at different redshifts
        ax1 = axes[0]
        for i, (k, Pk) in enumerate(self.power_spectra):
            a = self.scale_factors[i]
            z = 1.0/a - 1
            ax1.loglog(k[1:], Pk[1:], label=f'z={z:.1f}')
        
        ax1.set_xlabel('k [h/Mpc]')
        ax1.set_ylabel('P(k) [(Mpc/h)³]')
        ax1.set_title('HQIV Matter Power Spectrum')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Growth of structure
        ax2 = axes[1]
        k_idx = len(self.power_spectra[0][0]) // 3
        growth = [Pk[k_idx] for k, Pk in self.power_spectra]
        
        ax2.semilogy(self.scale_factors, growth, 'b-', linewidth=2)
        ax2.set_xlabel('Scale factor a')
        ax2.set_ylabel(f'P(k) at k≈{self.power_spectra[0][0][k_idx]:.2f} h/Mpc')
        ax2.set_title('Growth of Structure')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        print(f"Power spectrum plot saved to {output_file}")
        
        return fig


def main():
    """Run a test simulation."""
    sim = HQIVSimulator(
        N_particles=32**3,
        L_box=100.0,
        N_grid=32,
        a_start=0.1,
        a_end=1.0,
        n_steps=100
    )
    
    power_spectra, scale_factors = sim.run(
        save_snapshots=True,
        snapshot_interval=20
    )
    
    sim.plot_power_spectrum('hqiv_power_spectrum.png')
    
    print("\n=== HQIV PM Simulation Summary ===")
    print(f"Final scale factor: {scale_factors[-1]:.4f}")
    print(f"Final redshift: {1/scale_factors[-1] - 1:.2f}")
    
    if len(power_spectra) >= 2:
        P_initial = np.max(power_spectra[0][1])
        P_final = np.max(power_spectra[-1][1])
        growth_factor = np.sqrt(P_final / P_initial)
        a_growth = scale_factors[-1] / scale_factors[0]
        
        print(f"\nGrowth factor: {growth_factor:.2f}")
        print(f"Linear growth (a ratio): {a_growth:.2f}")
        print(f"Enhancement: {growth_factor / a_growth:.2f}x linear")
        
        # In ΛCDM, growth factor ~ a, so P(k) ~ a²
        # Growth factor should be ~ a_growth
        print(f"\nExpected ΛCDM growth: ~{a_growth:.2f}")
        print(f"HQIV shows: {growth_factor:.2f}")
    
    return sim


if __name__ == "__main__":
    sim = main()