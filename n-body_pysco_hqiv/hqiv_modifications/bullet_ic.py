"""
HQIV Bullet Cluster Initial Conditions Module
==============================================

Generates initial conditions for the Bullet Cluster collision test.

Paper Reference: paper/main.tex, Section 9
    Two NFW halos colliding at ~3000-4000 km/s
    Separate particle species: "gas" (collisional, high-a) and "galaxies" (collisionless, low-a)
    Output: X-ray-like gas map + synthetic weak-lensing convergence κ-map

The Bullet Cluster is a key test for modified gravity/inertia theories.
The observed separation between X-ray gas and weak-lensing mass peak
must be reproduced without dark matter.

Author: HQIV Team
"""

import numpy as np
import configparser
from pathlib import Path

# Physical constants (SI units)
c = 2.99792458e8       # Speed of light [m/s]
G0_SI = 6.67430e-11    # Gravitational constant [m³ kg⁻¹ s⁻²]
Mpc_m = 3.0856775814913673e22  # Megaparsec in meters
Msun = 1.98847e30      # Solar mass [kg]
kpc_m = 3.0856775814913673e19  # Kiloparsec in meters


def load_config():
    """Load configuration from config_hqiv.ini."""
    config_path = Path(__file__).parent / 'config_hqiv.ini'
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def nfw_density(r, rs, rho_s):
    """
    NFW density profile.
    
    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
    
    Parameters
    ----------
    r : float or array
        Radial distance [Mpc]
    rs : float
        Scale radius [Mpc]
    rho_s : float
        Characteristic density [Msun/Mpc³]
        
    Returns
    -------
    rho : float or array
        Density [Msun/Mpc³]
    """
    x = r / rs
    return rho_s / (x * (1 + x)**2)


def nfw_mass_enclosed(r, rs, rho_s):
    """
    Enclosed mass within radius r for NFW profile.
    
    M(<r) = 4π ρ_s r_s³ [ln(1 + r/r_s) - (r/r_s)/(1 + r/r_s)]
    
    Parameters
    ----------
    r : float or array
        Radial distance [Mpc]
    rs : float
        Scale radius [Mpc]
    rho_s : float
        Characteristic density [Msun/Mpc³]
        
    Returns
    -------
    M : float or array
        Enclosed mass [Msun]
    """
    x = r / rs
    return 4 * np.pi * rho_s * rs**3 * (np.log(1 + x) - x / (1 + x))


def nfw_concentration_to_rho_s(M200, c, rho_crit):
    """
    Compute NFW characteristic density from concentration and M200.
    
    Parameters
    ----------
    M200 : float
        Mass within R200 [Msun]
    c : float
        Concentration parameter c = R200 / rs
    rho_crit : float
        Critical density [Msun/Mpc³]
        
    Returns
    -------
    rho_s : float
        Characteristic density [Msun/Mpc³]
    """
    # Δ_c = 200 for M200 definition
    delta_c = 200 * c**3 / (3 * (np.log(1 + c) - c / (1 + c)))
    rho_s = delta_c * rho_crit
    return rho_s


def generate_nfw_halo(M200, c, N_particles, box_size, center, velocity,
                       r_max=None, seed=None):
    """
    Generate particles following an NFW profile.
    
    Parameters
    ----------
    M200 : float
        Mass within R200 [Msun]
    c : float
        Concentration parameter
    N_particles : int
        Number of particles
    box_size : float
        Box size [Mpc]
    center : array
        Center position [Mpc]
    velocity : array
        Bulk velocity [km/s]
    r_max : float, optional
        Maximum radius [Mpc]. If None, use R200.
    seed : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    positions : array
        Particle positions [N_particles, 3] [Mpc]
    velocities : array
        Particle velocities [N_particles, 3] [km/s]
    masses : array
        Particle masses [Msun]
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Critical density (assuming H0 = 73.2 km/s/Mpc)
    H0 = 73.2 * 1000.0 / Mpc_m  # 1/s
    rho_crit = 3 * H0**2 / (8 * np.pi * G0_SI) * Mpc_m**3 / Msun  # Msun/Mpc³
    
    # Scale radius
    R200 = (3 * M200 / (4 * np.pi * 200 * rho_crit))**(1/3)
    rs = R200 / c
    
    if r_max is None:
        r_max = R200
    
    # Characteristic density
    rho_s = nfw_concentration_to_rho_s(M200, c, rho_crit)
    
    # Generate radii using inverse transform sampling
    # CDF for NFW: M(<r) / M(<r_max)
    M_max = nfw_mass_enclosed(r_max, rs, rho_s)
    
    # Use rejection sampling for radii
    positions = np.zeros((N_particles, 3))
    
    # Generate random points in a sphere and accept/reject based on NFW density
    r_samples = np.random.uniform(0, r_max, N_particles * 10)
    theta_samples = np.random.uniform(0, np.pi, N_particles * 10)
    phi_samples = np.random.uniform(0, 2*np.pi, N_particles * 10)
    
    # Density at sample points
    rho_samples = nfw_density(r_samples, rs, rho_s)
    rho_max = rho_s / 0.25  # Maximum density at r = rs/2 approximately
    
    # Accept/reject
    u = np.random.uniform(0, 1, N_particles * 10)
    accept = u < rho_samples / rho_max
    
    # Take first N_particles accepted
    r_accepted = r_samples[accept][:N_particles]
    theta_accepted = theta_samples[accept][:N_particles]
    phi_accepted = phi_samples[accept][:N_particles]
    
    # If not enough, generate more
    while len(r_accepted) < N_particles:
        r_new = np.random.uniform(0, r_max, N_particles)
        theta_new = np.random.uniform(0, np.pi, N_particles)
        phi_new = np.random.uniform(0, 2*np.pi, N_particles)
        rho_new = nfw_density(r_new, rs, rho_s)
        u_new = np.random.uniform(0, 1, N_particles)
        accept_new = u_new < rho_new / rho_max
        r_accepted = np.concatenate([r_accepted, r_new[accept_new]])[:N_particles]
        theta_accepted = np.concatenate([theta_accepted, theta_new[accept_new]])[:N_particles]
        phi_accepted = np.concatenate([phi_accepted, phi_new[accept_new]])[:N_particles]
    
    # Convert to Cartesian
    positions[:, 0] = r_accepted * np.sin(theta_accepted) * np.cos(phi_accepted)
    positions[:, 1] = r_accepted * np.sin(theta_accepted) * np.sin(phi_accepted)
    positions[:, 2] = r_accepted * np.cos(theta_accepted)
    
    # Add center offset
    positions += center
    
    # Apply periodic boundary conditions
    positions = positions % box_size
    
    # Velocities: bulk velocity + velocity dispersion
    # Velocity dispersion from virial theorem
    sigma_v = np.sqrt(G0_SI * M200 * Msun / (R200 * Mpc_m)) / 1000  # km/s
    
    velocities = np.zeros((N_particles, 3))
    velocities[:, 0] = velocity[0] + np.random.normal(0, sigma_v * 0.3, N_particles)
    velocities[:, 1] = velocity[1] + np.random.normal(0, sigma_v * 0.3, N_particles)
    velocities[:, 2] = velocity[2] + np.random.normal(0, sigma_v * 0.3, N_particles)
    
    # Particle masses
    masses = np.full(N_particles, M200 / N_particles)
    
    return positions, velocities, masses


def generate_gas_particles(M_gas, N_particles, box_size, center, velocity,
                            T_gas=1e8, seed=None):
    """
    Generate gas particles with thermal distribution.
    
    Gas in galaxy clusters is hot (T ~ 10^8 K) and follows
    a hydrostatic equilibrium profile (simplified as isothermal beta-model).
    
    Parameters
    ----------
    M_gas : float
        Total gas mass [Msun]
    N_particles : int
        Number of particles
    box_size : float
        Box size [Mpc]
    center : array
        Center position [Mpc]
    velocity : array
        Bulk velocity [km/s]
    T_gas : float
        Gas temperature [K]
    seed : int, optional
        Random seed
        
    Returns
    -------
    positions : array
        Particle positions [N_particles, 3] [Mpc]
    velocities : array
        Particle velocities [N_particles, 3] [km/s]
    masses : array
        Particle masses [Msun]
    thermal_accel : array
        Thermal acceleration magnitude for each particle [m/s²]
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Beta-model core radius (typical for clusters)
    r_core = 0.2  # Mpc
    r_max = 2.0   # Mpc
    
    # Beta-model: ρ(r) ∝ (1 + (r/r_c)²)^(-3β/2)
    # Typical β ≈ 2/3
    beta = 0.67
    
    # Generate radii using rejection sampling
    r_samples = np.random.uniform(0, r_max, N_particles * 10)
    rho_samples = (1 + (r_samples / r_core)**2)**(-1.5 * beta)
    rho_max = 1.0
    
    u = np.random.uniform(0, 1, N_particles * 10)
    accept = u < rho_samples / rho_max
    
    r_accepted = r_samples[accept][:N_particles]
    
    while len(r_accepted) < N_particles:
        r_new = np.random.uniform(0, r_max, N_particles)
        rho_new = (1 + (r_new / r_core)**2)**(-1.5 * beta)
        u_new = np.random.uniform(0, 1, N_particles)
        accept_new = u_new < rho_new / rho_max
        r_accepted = np.concatenate([r_accepted, r_new[accept_new]])[:N_particles]
    
    # Angular positions
    theta = np.random.uniform(0, np.pi, N_particles)
    phi = np.random.uniform(0, 2*np.pi, N_particles)
    
    # Convert to Cartesian
    positions = np.zeros((N_particles, 3))
    positions[:, 0] = r_accepted * np.sin(theta) * np.cos(phi)
    positions[:, 1] = r_accepted * np.sin(theta) * np.sin(phi)
    positions[:, 2] = r_accepted * np.cos(theta)
    
    # Add center offset
    positions += center
    
    # Apply periodic boundary conditions
    positions = positions % box_size
    
    # Thermal velocity dispersion
    # σ_v = sqrt(k_B T / m_p) for protons
    k_B = 1.380649e-23  # J/K
    m_p = 1.6726219e-27  # kg
    sigma_thermal = np.sqrt(k_B * T_gas / m_p) / 1000  # km/s
    
    velocities = np.zeros((N_particles, 3))
    velocities[:, 0] = velocity[0] + np.random.normal(0, sigma_thermal, N_particles)
    velocities[:, 1] = velocity[1] + np.random.normal(0, sigma_thermal, N_particles)
    velocities[:, 2] = velocity[2] + np.random.normal(0, sigma_thermal, N_particles)
    
    # Particle masses
    masses = np.full(N_particles, M_gas / N_particles)
    
    # Thermal acceleration (for inertia factor calculation)
    # a_thermal ~ k_B T / (m_p * r_core) in physical units
    thermal_accel = k_B * T_gas / (m_p * r_core * Mpc_m) * np.ones(N_particles)
    
    return positions, velocities, masses, thermal_accel


class BulletClusterIC:
    """
    Generate initial conditions for the Bullet Cluster collision.
    
    Parameters
    ----------
    config : configparser.ConfigParser, optional
        Configuration object. If None, load from default file.
    """
    
    def __init__(self, config=None):
        if config is None:
            config = load_config()
        
        self.config = config
        
        # Parse parameters
        self.M_main = float(config.get('bullet_cluster', 'M_main', fallback='2.5e14'))
        self.M_sub = float(config.get('bullet_cluster', 'M_sub', fallback='1.5e14'))
        self.v_collision = float(config.get('bullet_cluster', 'v_collision', fallback='4500'))
        self.impact_param = float(config.get('bullet_cluster', 'impact_param', fallback='150'))
        self.gas_fraction = float(config.get('bullet_cluster', 'gas_fraction', fallback='0.15'))
        self.c_main = float(config.get('bullet_cluster', 'c_main', fallback='4.0'))
        self.c_sub = float(config.get('bullet_cluster', 'c_sub', fallback='3.5'))
        self.initial_separation = float(config.get('bullet_cluster', 'initial_separation', fallback='2.0'))
        
        self.box_size = float(config.get('simulation', 'boxlen', fallback='100.0'))
        self.npart_gas = int(config.get('simulation', 'npart_gas', fallback='134217728'))
        self.npart_galaxies = int(config.get('simulation', 'npart_galaxies', fallback='16777216'))
        
        # Derived parameters
        self.M_main_gas = self.M_main * self.gas_fraction
        self.M_main_gal = self.M_main * (1 - self.gas_fraction)
        self.M_sub_gas = self.M_sub * self.gas_fraction
        self.M_sub_gal = self.M_sub * (1 - self.gas_fraction)
        
    def generate(self, seed=42):
        """
        Generate full initial conditions for the collision.
        
        Parameters
        ----------
        seed : int
            Random seed for reproducibility
            
        Returns
        -------
        ic_dict : dict
            Dictionary containing all particle data:
            - 'pos_gas': gas particle positions
            - 'vel_gas': gas particle velocities
            - 'mass_gas': gas particle masses
            - 'pos_gal': galaxy particle positions
            - 'vel_gal': galaxy particle velocities
            - 'mass_gal': galaxy particle masses
            - 'species': species labels
            - 'thermal_accel': thermal acceleration for gas
        """
        print("Generating Bullet Cluster initial conditions...")
        print(f"  Main cluster: M = {self.M_main:.2e} Msun")
        print(f"  Subcluster: M = {self.M_sub:.2e} Msun")
        print(f"  Collision velocity: {self.v_collision:.0f} km/s")
        
        # Center of box
        center = np.array([self.box_size / 2, self.box_size / 2, self.box_size / 2])
        
        # Initial positions (along x-axis)
        center_main = center + np.array([self.initial_separation / 2, 0, 0])
        center_sub = center - np.array([self.initial_separation / 2, self.impact_param / 1000, 0])
        
        # Velocities (approaching each other)
        v_main = np.array([-self.v_collision / 2, 0, 0])
        v_sub = np.array([self.v_collision / 2, 0, 0])
        
        # Number of particles per component
        N_gas_main = int(self.npart_gas * self.M_main_gas / (self.M_main_gas + self.M_sub_gas))
        N_gas_sub = self.npart_gas - N_gas_main
        N_gal_main = int(self.npart_galaxies * self.M_main_gal / (self.M_main_gal + self.M_sub_gal))
        N_gal_sub = self.npart_galaxies - N_gal_main
        
        print(f"  Gas particles: {N_gas_main} (main) + {N_gas_sub} (sub)")
        print(f"  Galaxy particles: {N_gal_main} (main) + {N_gal_sub} (sub)")
        
        # Generate main cluster
        print("\n  Generating main cluster...")
        
        pos_gal_main, vel_gal_main, mass_gal_main = generate_nfw_halo(
            self.M_main_gal, self.c_main, N_gal_main,
            self.box_size, center_main, v_main, seed=seed
        )
        
        pos_gas_main, vel_gas_main, mass_gas_main, thermal_main = generate_gas_particles(
            self.M_main_gas, N_gas_main,
            self.box_size, center_main, v_main, seed=seed + 1
        )
        
        # Generate subcluster (bullet)
        print("  Generating subcluster (bullet)...")
        
        pos_gal_sub, vel_gal_sub, mass_gal_sub = generate_nfw_halo(
            self.M_sub_gal, self.c_sub, N_gal_sub,
            self.box_size, center_sub, v_sub, seed=seed + 2
        )
        
        pos_gas_sub, vel_gas_sub, mass_gas_sub, thermal_sub = generate_gas_particles(
            self.M_sub_gas, N_gas_sub,
            self.box_size, center_sub, v_sub, seed=seed + 3
        )
        
        # Combine all particles
        pos_gas = np.vstack([pos_gas_main, pos_gas_sub])
        vel_gas = np.vstack([vel_gas_main, vel_gas_sub])
        mass_gas = np.concatenate([mass_gas_main, mass_gas_sub])
        thermal_accel = np.concatenate([thermal_main, thermal_sub])
        
        pos_gal = np.vstack([pos_gal_main, pos_gal_sub])
        vel_gal = np.vstack([vel_gal_main, vel_gal_sub])
        mass_gal = np.concatenate([mass_gal_main, mass_gal_sub])
        
        # Species labels
        species_gas = np.full(self.npart_gas, 'gas', dtype='<U10')
        species_gal = np.full(self.npart_galaxies, 'galaxy', dtype='<U10')
        
        print(f"\n  Total particles: {self.npart_gas + self.npart_galaxies}")
        print(f"  Total mass: {self.M_main + self.M_sub:.2e} Msun")
        
        return {
            'pos_gas': pos_gas,
            'vel_gas': vel_gas,
            'mass_gas': mass_gas,
            'pos_gal': pos_gal,
            'vel_gal': vel_gal,
            'mass_gal': mass_gal,
            'species_gas': species_gas,
            'species_gal': species_gal,
            'thermal_accel': thermal_accel,
            'N_total': self.npart_gas + self.npart_galaxies,
            'M_total': self.M_main + self.M_sub,
        }
    
    def get_collision_parameters(self):
        """Return collision parameters for analysis."""
        return {
            'M_main': self.M_main,
            'M_sub': self.M_sub,
            'v_collision': self.v_collision,
            'impact_param': self.impact_param,
            'initial_separation': self.initial_separation,
        }


def generate_collision_ic(config=None, seed=42):
    """
    Convenience function to generate Bullet Cluster IC.
    
    Parameters
    ----------
    config : configparser.ConfigParser, optional
        Configuration object
    seed : int
        Random seed
        
    Returns
    -------
    ic_dict : dict
        Initial conditions dictionary
    """
    ic_gen = BulletClusterIC(config)
    return ic_gen.generate(seed=seed)


# =============================================================================
# Testing
# =============================================================================

def test_bullet_ic():
    """Test Bullet Cluster initial conditions generation."""
    import matplotlib.pyplot as plt
    
    print("Testing Bullet Cluster IC generation...")
    
    # Use smaller particle counts for testing
    test_config = configparser.ConfigParser()
    test_config.read_dict({
        'bullet_cluster': {
            'M_main': '2.5e14',
            'M_sub': '1.5e14',
            'v_collision': '4500',
            'impact_param': '150',
            'gas_fraction': '0.15',
            'c_main': '4.0',
            'c_sub': '3.5',
            'initial_separation': '2.0',
        },
        'simulation': {
            'boxlen': '50.0',
            'npart_gas': '10000',
            'npart_galaxies': '5000',
        }
    })
    
    ic = generate_collision_ic(test_config, seed=42)
    
    print(f"\nGenerated IC summary:")
    print(f"  Gas particles: {len(ic['pos_gas'])}")
    print(f"  Galaxy particles: {len(ic['pos_gal'])}")
    print(f"  Position ranges:")
    print(f"    Gas: [{ic['pos_gas'].min():.2f}, {ic['pos_gas'].max():.2f}] Mpc")
    print(f"    Gal: [{ic['pos_gal'].min():.2f}, {ic['pos_gal'].max():.2f}] Mpc")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    ax = axes[0]
    ax.scatter(ic['pos_gas'][:, 0], ic['pos_gas'][:, 1], s=0.1, alpha=0.5, label='Gas')
    ax.scatter(ic['pos_gal'][:, 0], ic['pos_gal'][:, 1], s=0.5, alpha=0.5, label='Galaxies', color='red')
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    ax.set_title('Bullet Cluster IC (x-y projection)')
    ax.legend()
    ax.set_aspect('equal')
    
    ax = axes[1]
    ax.scatter(ic['pos_gas'][:, 0], ic['pos_gas'][:, 2], s=0.1, alpha=0.5, label='Gas')
    ax.scatter(ic['pos_gal'][:, 0], ic['pos_gal'][:, 2], s=0.5, alpha=0.5, label='Galaxies', color='red')
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('z [Mpc]')
    ax.set_title('Bullet Cluster IC (x-z projection)')
    ax.legend()
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig('bullet_ic_test.png', dpi=150)
    plt.close()
    print("\nSaved bullet_ic_test.png")
    
    return ic


if __name__ == "__main__":
    test_bullet_ic()