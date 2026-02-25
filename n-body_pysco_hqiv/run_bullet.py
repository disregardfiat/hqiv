#!/usr/bin/env python3
"""
HQIV Bullet Cluster Simulation
==============================

Main script to run the Bullet Cluster N-body simulation with HQIV physics.

Paper Reference: paper/main.tex, Section 9
    Tests the HQIV framework against the Bullet Cluster observations:
    - Two NFW halos colliding at ~3000-4000 km/s
    - Separate particle species: gas (high-a) and galaxies (low-a)
    - Output: X-ray-like gas map + synthetic weak-lensing κ-map
    - Compare offset (~180 kpc) and lensing/gas ratio (~4-6:1) to real data

Usage:
    python run_bullet.py --resolution 256 --npart 1e7

Author: HQIV Team
"""

import argparse
import configparser
import os
import sys
import time
from pathlib import Path

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq

# Add hqiv_modifications to path
sys.path.insert(0, str(Path(__file__).parent))

# Numba parallel CIC (used if numba available)
try:
    from numba import njit, prange
    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

if _NUMBA_AVAILABLE:
    @njit(parallel=True, fastmath=True)
    def _cic_density_parallel(positions, density_tls, dx, N_grid, nthreads):
        """Thread-local CIC deposit; sum density_tls over axis 0 to get density."""
        npart = positions.shape[0]
        chunk = (npart + nthreads - 1) // nthreads
        for tid in prange(nthreads):
            start = tid * chunk
            end = min(start + chunk, npart)
            for i in range(start, end):
                x, y, z = positions[i, 0], positions[i, 1], positions[i, 2]
                ix = int(np.floor(x / dx)) % N_grid
                iy = int(np.floor(y / dx)) % N_grid
                iz = int(np.floor(z / dx)) % N_grid
                ux = x / dx - np.floor(x / dx)
                uy = y / dx - np.floor(y / dx)
                uz = z / dx - np.floor(z / dx)
                for dx_ in (0, 1):
                    for dy in (0, 1):
                        for dz in (0, 1):
                            w = ((1 - ux) if dx_ == 0 else ux) * (
                                (1 - uy) if dy == 0 else uy
                            ) * ((1 - uz) if dz == 0 else uz)
                            ixc = (ix + dx_) % N_grid
                            iyc = (iy + dy) % N_grid
                            izc = (iz + dz) % N_grid
                            density_tls[tid, ixc, iyc, izc] += w

    @njit(parallel=True, fastmath=True)
    def _cic_velocity_parallel(positions, velocities, vx_tls, vy_tls, vz_tls, w_tls, dx, N_grid, nthreads):
        """Thread-local CIC velocity; sum then normalize outside."""
        npart = positions.shape[0]
        chunk = (npart + nthreads - 1) // nthreads
        for tid in prange(nthreads):
            start = tid * chunk
            end = min(start + chunk, npart)
            for i in range(start, end):
                x, y, z = positions[i, 0], positions[i, 1], positions[i, 2]
                vx_p, vy_p, vz_p = velocities[i, 0], velocities[i, 1], velocities[i, 2]
                ix = int(np.floor(x / dx)) % N_grid
                iy = int(np.floor(y / dx)) % N_grid
                iz = int(np.floor(z / dx)) % N_grid
                ux = x / dx - np.floor(x / dx)
                uy = y / dx - np.floor(y / dx)
                uz = z / dx - np.floor(z / dx)
                for dx_ in (0, 1):
                    for dy in (0, 1):
                        for dz in (0, 1):
                            w = ((1 - ux) if dx_ == 0 else ux) * (
                                (1 - uy) if dy == 0 else uy
                            ) * ((1 - uz) if dz == 0 else uz)
                            ixc, iyc, izc = (ix + dx_) % N_grid, (iy + dy) % N_grid, (iz + dz) % N_grid
                            vx_tls[tid, ixc, iyc, izc] += w * vx_p
                            vy_tls[tid, ixc, iyc, izc] += w * vy_p
                            vz_tls[tid, ixc, iyc, izc] += w * vz_p
                            w_tls[tid, ixc, iyc, izc] += w

    @njit(parallel=True, fastmath=True)
    def _interpolate_to_particles_parallel(positions, field, dx, N_grid):
        """CIC interpolate field to particle positions. field shape (3, N, N, N)."""
        npart = positions.shape[0]
        values = np.empty((npart, 3))
        for i in prange(npart):
            x, y, z = positions[i, 0], positions[i, 1], positions[i, 2]
            ix = int(np.floor(x / dx)) % N_grid
            iy = int(np.floor(y / dx)) % N_grid
            iz = int(np.floor(z / dx)) % N_grid
            ux = x / dx - np.floor(x / dx)
            uy = y / dx - np.floor(y / dx)
            uz = z / dx - np.floor(z / dx)
            v0 = v1 = v2 = 0.0
            for dx_ in (0, 1):
                for dy in (0, 1):
                    for dz in (0, 1):
                        w = ((1 - ux) if dx_ == 0 else ux) * (
                            (1 - uy) if dy == 0 else uy
                        ) * ((1 - uz) if dz == 0 else uz)
                        ixc, iyc, izc = (ix + dx_) % N_grid, (iy + dy) % N_grid, (iz + dz) % N_grid
                        v0 += w * field[0, ixc, iyc, izc]
                        v1 += w * field[1, ixc, iyc, izc]
                        v2 += w * field[2, ixc, iyc, izc]
            values[i, 0], values[i, 1], values[i, 2] = v0, v1, v2
        return values
else:
    _cic_density_parallel = _cic_velocity_parallel = _interpolate_to_particles_parallel = None

from hqiv_modifications import (
    compute_phi_field,
    inertia_reduction_factor,
    vorticity_source_term,
    effective_gravitational_constant,
    BulletClusterIC,
    generate_collision_ic,
)
from hqiv_modifications.phi_field import PhiFieldSolver
from hqiv_modifications.inertia_factor import InertiaFactorSolver
from hqiv_modifications.vorticity_source import VorticitySolver
from hqiv_modifications.g_eff import EffectiveGravitySolver

# Physical constants
c = 2.99792458e8       # Speed of light [m/s]
G0_SI = 6.67430e-11    # Gravitational constant [m³ kg⁻¹ s⁻²]
Mpc_m = 3.0856775814913673e22  # Megaparsec in meters
Msun = 1.98847e30      # Solar mass [kg]


class HQIVBulletSimulation:
    """
    Full HQIV N-body simulation for the Bullet Cluster.
    
    This class orchestrates:
    1. Initial condition generation
    2. PM force calculation with HQIV modifications
    3. Time integration with inertia reduction
    4. Vorticity evolution
    5. Output generation
    
    Parameters
    ----------
    config_path : str
        Path to configuration file
    resolution : int
        Grid resolution (power of 2)
    npart : int
        Approximate total number of particles
    output_dir : str
        Output directory
    """
    
    def __init__(self, config_path=None, resolution=256, npart=1e7, output_dir='./output/',
                 gamma=None, omegab=None, omegam=None, h0=None, alpha=None, z=None, box_size=None,
                 nthreads=None):
        # Load configuration
        if config_path is None:
            config_path = Path(__file__).parent / 'hqiv_modifications' / 'config_hqiv.ini'
        
        self.config = configparser.ConfigParser()
        self.config.read(config_path)
        
        # Override with command-line arguments
        self.N_grid = resolution
        self.npart_total = int(npart)
        self.nthreads = int(nthreads) if nthreads is not None else int(self.config.get('simulation', 'nthreads', fallback='4'))
        # Thread-local CIC uses nthreads_cic to limit memory (each thread has N^3 arrays)
        self._nthreads_cic = min(self.nthreads, 8) if _NUMBA_AVAILABLE else 1
        
        # Parse configuration
        self.H0 = float(self.config.get('cosmology', 'H0', fallback='73.2'))
        self.Omega_m = float(self.config.get('cosmology', 'Omega_m', fallback='0.06'))
        self.gamma = float(self.config.get('hqiv_parameters', 'gamma', fallback='0.40'))
        self.alpha_G = float(self.config.get('hqiv_parameters', 'alpha_G', fallback='0.6'))
        self.chi = float(self.config.get('hqiv_parameters', 'chi', fallback='0.172'))
        self.f_min = float(self.config.get('hqiv_parameters', 'f_min', fallback='0.01'))
        
        self.box_size = float(self.config.get('simulation', 'boxlen', fallback='100.0'))
        if box_size is not None:
            self.box_size = float(box_size)
        self.a_start = float(self.config.get('simulation', 'a_start', fallback='0.5'))
        self.a_end = float(self.config.get('simulation', 'a_end', fallback='1.0'))
        self.n_steps = int(self.config.get('simulation', 'n_steps', fallback='500'))
        
        # CLI overrides (e.g. peak-alignment best-fit)
        if omegam is not None:
            self.Omega_m = float(omegam)
        elif omegab is not None:
            self.Omega_m = float(omegab)
        if h0 is not None:
            self.H0 = float(h0) * 100.0  # dimensionless h -> H0 in km/s/Mpc
        if gamma is not None:
            self.gamma = float(gamma)
        if alpha is not None:
            self.alpha_G = float(alpha)
        if z is not None:
            self.a_end = 1.0 / (1.0 + float(z))
        
        # Output
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Derived quantities
        self.dx = self.box_size / self.N_grid
        
        # Initialize solvers
        self._init_solvers()
        
        # Storage
        self.positions = None
        self.velocities = None
        self.masses = None
        self.species = None
        self.thermal_accel = None
        
        self.time = 0.0
        self.a_current = self.a_start
        
        # Diagnostics
        self.history = {
            'a': [],
            'time': [],
            'offset_gas_gal': [],
            'offset_lensing_gas': [],
            'vorticity_rms': [],
            'G_eff': [],
        }
        
    def _init_solvers(self):
        """Initialize HQIV physics solvers."""
        print("Initializing HQIV solvers...")
        
        # Background cosmology
        self.H0_SI = self.H0 * 1000.0 / Mpc_m
        
        # Phi field solver
        self.phi_solver = PhiFieldSolver(
            H0=self.H0, Omega_m=self.Omega_m, gamma=self.gamma,
            N_grid=self.N_grid, box_size=self.box_size
        )
        
        # Inertia factor solver
        self.inertia_solver = InertiaFactorSolver(
            H0=self.H0, chi=self.chi, f_min=self.f_min
        )
        
        # Vorticity solver
        self.vorticity_solver = VorticitySolver(
            H0=self.H0, chi=self.chi, f_min=self.f_min,
            N_grid=self.N_grid, box_size=self.box_size
        )
        
        # Effective gravity solver
        self.gravity_solver = EffectiveGravitySolver(
            H0=self.H0, Omega_m=self.Omega_m, gamma=self.gamma,
            alpha_G=self.alpha_G, N_grid=self.N_grid, box_size=self.box_size
        )
        
        print(f"  Grid: {self.N_grid}³ = {self.N_grid**3:,} cells")
        print(f"  Box: {self.box_size} Mpc")
        print(f"  dx = {self.dx:.4f} Mpc")
        
    def generate_ic(self, seed=42):
        """Generate initial conditions for the Bullet Cluster."""
        print("\n" + "="*60)
        print("Generating Initial Conditions")
        print("="*60)
        
        # Adjust particle counts based on npart_total
        gas_fraction = 0.8  # 80% gas particles
        npart_gas = int(self.npart_total * gas_fraction)
        npart_gal = self.npart_total - npart_gas
        
        # Update config
        self.config.set('simulation', 'npart_gas', str(npart_gas))
        self.config.set('simulation', 'npart_galaxies', str(npart_gal))
        self.config.set('simulation', 'boxlen', str(self.box_size))
        
        # Generate IC
        ic = generate_collision_ic(self.config, seed=seed)
        
        # Combine gas and galaxy particles
        self.positions = np.vstack([ic['pos_gas'], ic['pos_gal']])
        self.velocities = np.vstack([ic['vel_gas'], ic['vel_gal']])
        self.masses = np.concatenate([ic['mass_gas'], ic['mass_gal']])
        self.species = np.concatenate([ic['species_gas'], ic['species_gal']])
        self.thermal_accel = np.concatenate([
            ic['thermal_accel'],
            np.zeros(len(ic['pos_gal']))  # No thermal acceleration for galaxies
        ])
        
        self.npart_gas = len(ic['pos_gas'])
        self.npart_gal = len(ic['pos_gal'])
        
        print(f"\n  Total particles: {len(self.positions):,}")
        print(f"    Gas: {self.npart_gas:,}")
        print(f"    Galaxies: {self.npart_gal:,}")
        
        return ic
    
    def compute_density(self):
        """Compute density field on grid using CIC (parallel when numba available)."""
        if _NUMBA_AVAILABLE and _cic_density_parallel is not None and self._nthreads_cic > 1:
            density_tls = np.zeros((self._nthreads_cic, self.N_grid, self.N_grid, self.N_grid), dtype=np.float64)
            _cic_density_parallel(
                self.positions, density_tls, self.dx, self.N_grid, self._nthreads_cic
            )
            density = density_tls.sum(axis=0)
        else:
            density = np.zeros((self.N_grid, self.N_grid, self.N_grid))
            for i in range(len(self.positions)):
                x, y, z = self.positions[i]
                ix = int(np.floor(x / self.dx)) % self.N_grid
                iy = int(np.floor(y / self.dx)) % self.N_grid
                iz = int(np.floor(z / self.dx)) % self.N_grid
                ux = x / self.dx - np.floor(x / self.dx)
                uy = y / self.dx - np.floor(y / self.dx)
                uz = z / self.dx - np.floor(z / self.dx)
                for dx in [0, 1]:
                    for dy in [0, 1]:
                        for dz in [0, 1]:
                            w = ((1-ux) if dx==0 else ux) * ((1-uy) if dy==0 else uy) * ((1-uz) if dz==0 else uz)
                            density[(ix+dx)%self.N_grid, (iy+dy)%self.N_grid, (iz+dz)%self.N_grid] += w
        rho_mean = len(self.positions) / self.N_grid**3
        delta = density / rho_mean - 1.0
        return delta, density
    
    def compute_velocity_field(self):
        """Compute velocity field on grid (parallel when numba available)."""
        if _NUMBA_AVAILABLE and _cic_velocity_parallel is not None and self._nthreads_cic > 1:
            nt = self._nthreads_cic
            vx_tls = np.zeros((nt, self.N_grid, self.N_grid, self.N_grid), dtype=np.float64)
            vy_tls = np.zeros_like(vx_tls)
            vz_tls = np.zeros_like(vx_tls)
            w_tls = np.zeros_like(vx_tls)
            _cic_velocity_parallel(
                self.positions, self.velocities,
                vx_tls, vy_tls, vz_tls, w_tls, self.dx, self.N_grid, nt
            )
            vx = vx_tls.sum(axis=0)
            vy = vy_tls.sum(axis=0)
            vz = vz_tls.sum(axis=0)
            weights = w_tls.sum(axis=0)
        else:
            vx = np.zeros((self.N_grid, self.N_grid, self.N_grid))
            vy = np.zeros_like(vx)
            vz = np.zeros_like(vx)
            weights = np.zeros_like(vx)
            for i in range(len(self.positions)):
                x, y, z = self.positions[i]
                vx_p, vy_p, vz_p = self.velocities[i]
                ix = int(np.floor(x / self.dx)) % self.N_grid
                iy = int(np.floor(y / self.dx)) % self.N_grid
                iz = int(np.floor(z / self.dx)) % self.N_grid
                ux = x / self.dx - np.floor(x / self.dx)
                uy = y / self.dx - np.floor(y / self.dx)
                uz = z / self.dx - np.floor(z / self.dx)
                for ddx in [0, 1]:
                    for ddy in [0, 1]:
                        for ddz in [0, 1]:
                            w = ((1-ux) if ddx==0 else ux) * ((1-uy) if ddy==0 else uy) * ((1-uz) if ddz==0 else uz)
                            vx[(ix+ddx)%self.N_grid, (iy+ddy)%self.N_grid, (iz+ddz)%self.N_grid] += w * vx_p
                            vy[(ix+ddx)%self.N_grid, (iy+ddy)%self.N_grid, (iz+ddz)%self.N_grid] += w * vy_p
                            vz[(ix+ddx)%self.N_grid, (iy+ddy)%self.N_grid, (iz+ddz)%self.N_grid] += w * vz_p
                            weights[(ix+ddx)%self.N_grid, (iy+ddy)%self.N_grid, (iz+ddz)%self.N_grid] += w
        mask = weights > 0
        vx[mask] /= weights[mask]
        vy[mask] /= weights[mask]
        vz[mask] /= weights[mask]
        return np.stack([vx, vy, vz])
    
    def step(self, a, da):
        """
        Advance simulation by one time step.
        
        Implements the HQIV-modified equations of motion:
        1. Compute density and velocity fields
        2. Compute φ field from expansion scalar
        3. Solve modified Poisson equation with G_eff and horizon term
        4. Apply inertia reduction to acceleration
        5. Update velocities and positions
        6. Apply vorticity amplification
        """
        H_a = self.phi_solver.H_of_a(a)
        dt = da / (a * H_a)
        
        # Compute density field
        delta, density = self.compute_density()
        
        # Compute velocity field
        velocity_field = self.compute_velocity_field()
        
        # Update φ field
        phi = self.phi_solver.update(delta, velocity_field, a)
        phi_particles = self.phi_solver.get_phi_at_positions(self.positions)
        
        # Solve for gravitational acceleration with HQIV modifications
        _, acceleration_field, G_eff = self.gravity_solver.full_poisson_solve(
            delta, a, H_a, include_horizon=True
        )
        
        # Interpolate acceleration to particles
        acc = self._interpolate_to_particles(acceleration_field)
        
        # Compute local acceleration magnitude
        acc_mag = np.sqrt(np.sum(acc**2, axis=1))
        
        # Convert to physical units for inertia factor
        # acc is in code units [H0 * Mpc], need [m/s²]
        # a_physical = a_code * H0 * Mpc_m / a (peculiar acceleration)
        acc_physical = acc_mag * self.H0_SI * Mpc_m / a
        
        # Compute inertia factor for each particle
        f_inertia = self.inertia_solver.compute_for_particles(
            acc_physical, phi_particles, self.species
        )
        
        # Apply inertia modification to acceleration
        acc_modified = self.inertia_solver.apply_to_acceleration(acc, f_inertia)
        
        # Limit acceleration for stability
        acc_max = 10.0
        acc_mag_mod = np.sqrt(np.sum(acc_modified**2, axis=1))
        too_large = acc_mag_mod > acc_max
        if np.any(too_large):
            acc_modified[too_large] *= acc_max / acc_mag_mod[too_large, np.newaxis]
        
        # Update velocities (kick-drift-kick leapfrog)
        # dv/da = -v/a + a_acc / (a * H)
        H_prime = H_a / (self.H0 * 1000.0 / Mpc_m)  # Dimensionless
        
        # Half kick
        self.velocities += 0.5 * da * (-self.velocities / a + acc_modified / (a * H_prime))
        
        # Drift
        self.positions += da * self.velocities / (H_prime * a * a)
        self.positions = self.positions % self.box_size
        
        # Recompute at new position
        a_new = a + da
        H_new = self.phi_solver.H_of_a(a_new)
        H_prime_new = H_new / (self.H0 * 1000.0 / Mpc_m)
        
        delta_new, _ = self.compute_density()
        velocity_field_new = self.compute_velocity_field()
        phi_new = self.phi_solver.update(delta_new, velocity_field_new, a_new)
        
        _, acc_field_new, _ = self.gravity_solver.full_poisson_solve(
            delta_new, a_new, H_new, include_horizon=True
        )
        acc_new = self._interpolate_to_particles(acc_field_new)
        
        acc_mag_new = np.sqrt(np.sum(acc_new**2, axis=1))
        acc_physical_new = acc_mag_new * self.H0_SI * Mpc_m / a_new
        phi_particles_new = self.phi_solver.get_phi_at_positions(self.positions)
        
        f_inertia_new = self.inertia_solver.compute_for_particles(
            acc_physical_new, phi_particles_new, self.species
        )
        acc_modified_new = self.inertia_solver.apply_to_acceleration(acc_new, f_inertia_new)
        
        acc_mag_mod_new = np.sqrt(np.sum(acc_modified_new**2, axis=1))
        too_large_new = acc_mag_mod_new > acc_max
        if np.any(too_large_new):
            acc_modified_new[too_large_new] *= acc_max / acc_mag_mod_new[too_large_new, np.newaxis]
        
        # Second half kick
        self.velocities += 0.5 * da * (-self.velocities / a_new + acc_modified_new / (a_new * H_prime_new))
        
        # Limit velocities
        vel_max = 100.0
        vel_mag = np.sqrt(np.sum(self.velocities**2, axis=1))
        too_fast = vel_mag > vel_max
        if np.any(too_fast):
            self.velocities[too_fast] *= vel_max / vel_mag[too_fast, np.newaxis]
        
        # Update time
        self.time += dt
        self.a_current = a_new
        
        return G_eff, np.mean(f_inertia)
    
    def _interpolate_to_particles(self, field):
        """Interpolate field to particle positions using CIC (parallel when numba available)."""
        if _NUMBA_AVAILABLE and _interpolate_to_particles_parallel is not None:
            return _interpolate_to_particles_parallel(self.positions, field, self.dx, self.N_grid)
        values = np.zeros((len(self.positions), 3))
        for i, pos in enumerate(self.positions):
            ix = int(np.floor(pos[0] / self.dx)) % self.N_grid
            iy = int(np.floor(pos[1] / self.dx)) % self.N_grid
            iz = int(np.floor(pos[2] / self.dx)) % self.N_grid
            ux = pos[0] / self.dx - np.floor(pos[0] / self.dx)
            uy = pos[1] / self.dx - np.floor(pos[1] / self.dx)
            uz = pos[2] / self.dx - np.floor(pos[2] / self.dx)
            for ddx in [0, 1]:
                for ddy in [0, 1]:
                    for ddz in [0, 1]:
                        w = ((1-ux) if ddx==0 else ux) * ((1-uy) if ddy==0 else uy) * ((1-uz) if ddz==0 else uz)
                        for d in range(3):
                            values[i, d] += w * field[d, (ix+ddx)%self.N_grid, (iy+ddy)%self.N_grid, (iz+ddz)%self.N_grid]
        return values
    
    def compute_diagnostics(self):
        """Compute diagnostic quantities."""
        # Separate gas and galaxy positions
        pos_gas = self.positions[:self.npart_gas]
        pos_gal = self.positions[self.npart_gas:]
        
        # Centers of mass
        com_gas = np.mean(pos_gas, axis=0)
        com_gal = np.mean(pos_gal, axis=0)
        
        # Offset between gas and galaxy centers
        offset_gas_gal = np.sqrt(np.sum((com_gas - com_gal)**2)) * 1000  # kpc
        
        return {
            'offset_gas_gal': offset_gas_gal,
            'com_gas': com_gas,
            'com_gal': com_gal,
        }
    
    def run(self, save_interval=10, verbose=True):
        """
        Run the full simulation.
        
        Parameters
        ----------
        save_interval : int
            Save snapshots every N steps
        verbose : bool
            Print progress
        """
        print("\n" + "="*60)
        print("Running HQIV Bullet Cluster Simulation")
        print("="*60)
        print(f"  a_start = {self.a_start:.4f}")
        print(f"  a_end = {self.a_end:.4f}")
        print(f"  n_steps = {self.n_steps}")
        print(f"  Output: {self.output_dir}")
        print("="*60)
        
        start_time = time.time()
        
        a_arr = np.linspace(self.a_start, self.a_end, self.n_steps + 1)
        
        for i in range(self.n_steps):
            a = a_arr[i]
            da = a_arr[i+1] - a_arr[i]
            
            G_eff, f_mean = self.step(a, da)
            
            # Record history
            self.history['a'].append(self.a_current)
            self.history['time'].append(self.time)
            self.history['G_eff'].append(G_eff)
            
            # Compute diagnostics
            if (i + 1) % save_interval == 0:
                diag = self.compute_diagnostics()
                self.history['offset_gas_gal'].append(diag['offset_gas_gal'])
                
                if verbose:
                    elapsed = time.time() - start_time
                    print(f"  Step {i+1}/{self.n_steps}: a={self.a_current:.4f}, "
                          f"z={1/self.a_current-1:.2f}, "
                          f"offset={diag['offset_gas_gal']:.1f} kpc, "
                          f"G_eff={G_eff:.3f}, "
                          f"<f>={f_mean:.3f}, "
                          f"t={elapsed:.1f}s")
                
                # Save snapshot
                self.save_snapshot(i + 1)
        
        print("\n" + "="*60)
        print("Simulation Complete!")
        print("="*60)
        # Always compute final diagnostics regardless of save_interval
        diag = self.compute_diagnostics()
        self.history['offset_gas_gal'].append(diag['offset_gas_gal'])
        self.history['offset_lensing_gas'].append(diag.get('offset_lensing_gas', 0.0))
        self.history['vorticity_rms'].append(diag.get('vorticity_rms', 0.0))

        elapsed = time.time() - start_time
        print(f"  Total time: {elapsed:.1f} seconds")
        print(f"  Final offset (gas-gal): {self.history['offset_gas_gal'][-1]:.1f} kpc")
        
        return self.history
    
    def save_snapshot(self, step):
        """Save simulation snapshot."""
        snapshot = {
            'positions': self.positions,
            'velocities': self.velocities,
            'masses': self.masses,
            'species': self.species,
            'a': self.a_current,
            'time': self.time,
        }
        
        np.savez(
            self.output_dir / f'snapshot_{step:05d}.npz',
            **snapshot
        )
    
    def save_final_results(self):
        """Save final results and analysis."""
        # Save particle data
        np.savez(
            self.output_dir / 'final_particles.npz',
            positions=self.positions,
            velocities=self.velocities,
            masses=self.masses,
            species=self.species,
        )
        
        # Save history
        np.savez(
            self.output_dir / 'simulation_history.npz',
            **self.history
        )
        
        # Compute and save density maps
        delta, density = self.compute_density()
        
        np.savez(
            self.output_dir / 'density_maps.npz',
            delta=delta,
            density=density,
        )
        
        print(f"\nSaved results to {self.output_dir}")


def main():
    """Main entry point."""
    print("run_bullet.py starting...", flush=True)
    parser = argparse.ArgumentParser(description='HQIV Bullet Cluster Simulation')
    parser.add_argument('--resolution', type=int, default=256,
                        help='Grid resolution (power of 2)')
    parser.add_argument('--npart', type=float, default=1e7,
                        help='Approximate number of particles')
    parser.add_argument('--config', type=str, default=None,
                        help='Configuration file path')
    parser.add_argument('--output', type=str, default='./output/',
                        help='Output directory')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--steps', type=int, default=None,
                        help='Number of steps (overrides config)')
    # Cosmology / HQIV overrides (e.g. peak-alignment best-fit)
    parser.add_argument('--z', type=float, default=None,
                        help='Final redshift (sets a_end = 1/(1+z))')
    parser.add_argument('--gamma', type=float, default=None,
                        help='HQIV gamma (horizon term)')
    parser.add_argument('--omegab', type=float, default=None,
                        help='Baryon density Omega_b (used as Omega_m if --omegam not set)')
    parser.add_argument('--omegam', type=float, default=None,
                        help='Matter density Omega_m (overrides config and --omegab)')
    parser.add_argument('--h0', type=float, default=None,
                        help='Dimensionless Hubble h = H0/100 (e.g. 0.732)')
    parser.add_argument('--alpha', type=float, default=None,
                        help='HQIV alpha_G for G_eff(a) = (H/H0)^alpha')
    parser.add_argument('--box', type=float, default=None,
                        help='Box size [Mpc] (cluster-scale e.g. 5)')
    parser.add_argument('--nthreads', type=int, default=None,
                        help='Number of threads (default: all cores / config nthreads)')
    
    args = parser.parse_args()
    
    # Use all available cores by default for speed
    nthreads = args.nthreads
    if nthreads is None:
        try:
            nthreads = len(os.sched_getaffinity(0))
        except AttributeError:
            nthreads = os.cpu_count() or 4
    if _NUMBA_AVAILABLE:
        import numba
        numba.set_num_threads(nthreads)
        print(f"Using {nthreads} threads (numba)")
    
    # Create simulation
    sim = HQIVBulletSimulation(
        config_path=args.config,
        resolution=args.resolution,
        npart=args.npart,
        output_dir=args.output,
        gamma=args.gamma,
        omegab=args.omegab,
        omegam=args.omegam,
        h0=args.h0,
        alpha=args.alpha,
        z=args.z,
        box_size=args.box,
        nthreads=nthreads,
    )
    
    # Override steps if specified
    if args.steps is not None:
        sim.n_steps = args.steps
    
    # Print input parameters table
    print("\n" + "="*80)
    print("HQIV BULLET CLUSTER SIMULATION - INPUT PARAMETERS")
    print("="*80)
    print(f"{'Parameter':<25} {'Value':<15} {'Unit':<10} {'Source':<20}")
    print("-" * 80)
    print(f"{'Resolution':<25} {sim.N_grid:<15} {'cells':<10} {'Command line':<20}")
    print(f"{'Total Particles':<25} {sim.npart_total:,<15} {'':<10} {'Command line':<20}")
    print(f"{'Box Size':<25} {sim.box_size:<15.1f} {'Mpc':<10} {'Command line' if args.box is not None else 'Config file':<20}")
    print(f"{'Grid Spacing':<25} {sim.dx:<15.4f} {'Mpc':<10} {'Calculated':<20}")
    print(f"{'H0':<25} {sim.H0:<15.1f} {'km/s/Mpc':<10} {'Config file':<20}")
    print(f"{'Omega_m':<25} {sim.Omega_m:<15.3f} {'':<10} {'Config file':<20}")
    print(f"{'Gamma':<25} {sim.gamma:<15.3f} {'':<10} {'Config file':<20}")
    print(f"{'Alpha_G':<25} {sim.alpha_G:<15.3f} {'':<10} {'Config file':<20}")
    print(f"{'a_start':<25} {sim.a_start:<15.4f} {'':<10} {'Config file':<20}")
    print(f"{'a_end':<25} {sim.a_end:<15.4f} {'':<10} {'Config file':<20}")
    print(f"{'n_steps':<25} {sim.n_steps:<15} {'':<10} {'Config file':<20}")
    print(f"{'Output Dir':<25} {str(sim.output_dir):<15} {'':<10} {'Command line':<20}")
    print(f"{'nthreads':<25} {sim.nthreads:<15} {'':<10} {'Command line':<20}")
    print(f"{'Random Seed':<25} {args.seed:<15} {'':<10} {'Command line':<20}")
    print("="*80)
    
    # Generate initial conditions
    sim.generate_ic(seed=args.seed)
    
    # Run simulation
    history = sim.run(save_interval=10, verbose=True)
    
    # Save final results
    sim.save_final_results()
    
    return sim, history


if __name__ == "__main__":
    # Filter empty args (e.g. from shell "$@" or nohup) so argparse doesn't fail
    sys.argv = [a for a in sys.argv if isinstance(a, str) and a.strip() != ""]
    try:
        sys.stdout.flush()
        sys.stderr.flush()
        sim, history = main()
    except Exception as e:
        import traceback
        print("run_bullet.py failed:", e, flush=True)
        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        raise