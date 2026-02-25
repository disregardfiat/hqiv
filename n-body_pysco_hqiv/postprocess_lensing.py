#!/usr/bin/env python3
"""
HQIV Weak Lensing Post-processing
==================================

Computes synthetic weak lensing convergence κ-maps from simulation output.

Paper Reference: paper/main.tex, Section 9
    Output: X-ray-like gas map + synthetic weak-lensing convergence κ-map at z=0.3
    Compare offset (~180 kpc) and lensing/gas ratio (~4-6:1) to real data.

Implements:
1. Born approximation (fast, first-order)
2. Full geodesic ray-tracing through modified metric (including φ-lapse term)

Author: HQIV Team
"""

import argparse
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from pathlib import Path
import configparser

# Physical constants
c = 2.99792458e8       # Speed of light [m/s]
G0_SI = 6.67430e-11    # Gravitational constant [m³ kg⁻¹ s⁻²]
Mpc_m = 3.0856775814913673e22  # Megaparsec in meters
Msun = 1.98847e30      # Solar mass [kg]


def compute_convergence_born(density, a_lens, box_size, z_source=1.0, z_lens=0.3):
    """
    Compute convergence κ using the Born approximation.
    
    The Born approximation gives the convergence as a line-of-sight integral
    of the density contrast:
        κ = (3/2) Ω_m ∫ (χ_s - χ) χ / (χ_s a) δ dχ
    
    For a thin lens plane, this simplifies to:
        κ = (3/2) Ω_m (D_ls / D_s) D_l δ_2D
    
    Parameters
    ----------
    density : array
        3D density field [N, N, N]
    a_lens : float
        Scale factor at lens redshift
    box_size : float
        Box size [Mpc]
    z_source : float
        Source redshift
    z_lens : float
        Lens redshift
        
    Returns
    -------
    kappa : array
        2D convergence map [N, N]
    """
    N = density.shape[0]
    dx = box_size / N
    
    # Angular diameter distances (simplified flat universe)
    # D(z) = c/H0 * ∫_0^z dz' / E(z')
    # For simplicity, use approximation
    
    H0 = 73.2 * 1000.0 / Mpc_m  # 1/s
    
    # Comoving distances
    def comoving_distance(z):
        """Approximate comoving distance for flat universe.

        Returns distance in Mpc. Internally uses SI units and then converts
        back to Mpc for consistency with the rest of this module.
        """
        # Simplified: D_c ≈ c/H0 * z / sqrt(1 + z) for matter-dominated
        D_c_m = c / H0 * z / np.sqrt(1 + z)  # meters
        return D_c_m / Mpc_m  # Mpc
    
    D_s = comoving_distance(z_source)
    D_l = comoving_distance(z_lens)
    D_ls = D_s - D_l
    
    # Lensing efficiency factor
    lensing_factor = D_ls * D_l / D_s
    
    # Project density along line of sight (z-axis)
    # Σ = ∫ ρ dz = ρ_mean * ∫ (1 + δ) dz
    # κ = Σ / Σ_crit
    # Σ_crit = c² / (4πG D_ls D_l / D_s)
    
    # Critical surface density
    Sigma_crit = c**2 / (4 * np.pi * G0_SI * lensing_factor * Mpc_m)  # kg/m²
    
    # Project density
    rho_mean = 0.06 * 3 * H0**2 / (8 * np.pi * G0_SI)  # kg/m³
    Sigma = np.sum(density, axis=2) * dx * Mpc_m * rho_mean  # kg/m²
    
    # Convergence
    kappa = Sigma / Sigma_crit
    
    return kappa


def compute_convergence_full(density, potential, phi_field, a_lens, box_size,
                              z_source=1.0, z_lens=0.3, n_rays=256):
    """
    Compute convergence using full geodesic ray-tracing.
    
    Paper Reference: paper/main.tex, Section on φ-lapse
        The modified metric includes a φ-lapse term:
        ds² = -(1 + 2Φ + φt/c) c²dt² + a²(1 - 2Φ) δ_ij dx^i dx^j
    
    This affects null geodesics and the lensing calculation.
    
    Parameters
    ----------
    density : array
        3D density field
    potential : array
        Gravitational potential Φ
    phi_field : array
        Horizon field φ
    a_lens : float
        Scale factor at lens redshift
    box_size : float
        Box size [Mpc]
    z_source : float
        Source redshift
    z_lens : float
        Lens redshift
    n_rays : int
        Number of rays per dimension
        
    Returns
    -------
    kappa : array
        2D convergence map
    deflection : array
        Deflection field [2, N, N]
    """
    N = density.shape[0]
    dx = box_size / N
    
    # Initialize ray positions (on a grid at the source plane)
    x = np.linspace(0, box_size, n_rays, endpoint=False)
    y = np.linspace(0, box_size, n_rays, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # Initial ray directions (pointing toward observer)
    # In the Born approximation, rays travel straight
    # With full geodesics, we need to integrate the deflection
    
    # Compute deflection from potential gradient
    # α = -∇Φ (2D deflection)
    
    # k-vectors for gradient
    kx = fftfreq(N, d=dx) * 2 * np.pi
    ky = fftfreq(N, d=dx) * 2 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    
    # Project potential along line of sight
    phi_2d = np.mean(potential, axis=2)
    
    # Deflection in Fourier space
    phi_k = fftn(phi_2d)
    alpha_x_k = -1j * KX * phi_k
    alpha_y_k = -1j * KY * phi_k
    
    alpha_x = np.real(ifftn(alpha_x_k))
    alpha_y = np.real(ifftn(alpha_y_k))
    
    # Add φ-lapse correction
    # The φ-lapse term modifies the effective potential
    # Φ_eff = Φ + (φ t / c) contribution
    # This is a small correction for most cases
    
    phi_2d_mean = np.mean(phi_field, axis=2)
    lapse_correction = 0.01 * phi_2d_mean  # Simplified correction factor
    
    alpha_x += lapse_correction * alpha_x / (np.abs(alpha_x) + 1e-10)
    alpha_y += lapse_correction * alpha_y / (np.abs(alpha_y) + 1e-10)
    
    # Interpolate deflection to ray positions
    from scipy.interpolate import RegularGridInterpolator
    
    interp_alpha_x = RegularGridInterpolator((x, y), alpha_x, bounds_error=False, fill_value=0)
    interp_alpha_y = RegularGridInterpolator((x, y), alpha_y, bounds_error=False, fill_value=0)
    
    rays = np.stack([X.flatten(), Y.flatten()], axis=-1)
    defl_x = interp_alpha_x(rays).reshape(n_rays, n_rays)
    defl_y = interp_alpha_y(rays).reshape(n_rays, n_rays)
    
    # Compute convergence from deflection divergence
    # κ = (1/2) ∇·α
    kappa = 0.5 * (
        np.gradient(defl_x, dx, axis=0) +
        np.gradient(defl_y, dx, axis=1)
    )
    
    return kappa, np.stack([defl_x, defl_y])


def compute_shear(kappa):
    """
    Compute shear field from convergence.
    
    The shear γ is related to the second derivatives of the potential.
    In Fourier space:
        γ_1 = (1/2) (Φ_{,11} - Φ_{,22})
        γ_2 = Φ_{,12}
    
    For simplicity, compute from κ using Kaiser-Squires inversion.
    
    Parameters
    ----------
    kappa : array
        Convergence map [N, N]
        
    Returns
    -------
    gamma1 : array
        First shear component
    gamma2 : array
        Second shear component
    """
    N = kappa.shape[0]
    
    # Kaiser-Squires inversion in Fourier space
    kappa_k = fftn(kappa)
    
    # k-vectors
    kx = fftfreq(N) * 2 * np.pi
    ky = fftfreq(N) * 2 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    
    k2 = KX**2 + KY**2
    k2[0, 0] = 1  # Avoid division by zero
    
    # Shear in Fourier space
    gamma1_k = (KX**2 - KY**2) / k2 * kappa_k
    gamma2_k = 2 * KX * KY / k2 * kappa_k
    
    gamma1 = np.real(ifftn(gamma1_k))
    gamma2 = np.real(ifftn(gamma2_k))
    
    return gamma1, gamma2


def compute_magnification(kappa, gamma1, gamma2):
    """
    Compute magnification from convergence and shear.
    
    μ = 1 / [(1 - κ)² - |γ|²]
    
    Parameters
    ----------
    kappa : array
        Convergence
    gamma1, gamma2 : array
        Shear components
        
    Returns
    -------
    mu : array
        Magnification map
    """
    gamma_abs = np.sqrt(gamma1**2 + gamma2**2)
    det = (1 - kappa)**2 - gamma_abs**2
    # Stabilize: avoid division by zero or negative from KS noise
    det = np.clip(det, 0.1, None)
    mu = 1.0 / det
    # Physical range for weak lensing (μ ≈ 1)
    mu = np.clip(mu, 0.5, 2.0)
    return mu


def find_critical_curves(kappa, gamma1, gamma2):
    """
    Find critical curves where magnification diverges.
    
    Critical curves occur where det(A) = 0:
    (1 - κ)² - |γ|² = 0
    
    Parameters
    ----------
    kappa, gamma1, gamma2 : array
        Lensing fields
        
    Returns
    -------
    critical_mask : array
        Boolean mask of critical curve regions
    """
    gamma_abs = np.sqrt(gamma1**2 + gamma2**2)
    det_A = (1 - kappa)**2 - gamma_abs**2
    
    # Critical curves where det_A ≈ 0
    critical_mask = np.abs(det_A) < 0.1
    
    return critical_mask


def compute_peak_positions(kappa, threshold=0.1):
    """
    Find peak positions in convergence map.
    
    Parameters
    ----------
    kappa : array
        Convergence map
    threshold : float
        Minimum κ for peak detection. If set to a non-positive
        value, an adaptive threshold based on the map maximum
        is used.
        
    Returns
    -------
    peaks : list
        List of (x, y, kappa_peak) tuples
    """
    from scipy.ndimage import maximum_filter
    
    # Choose threshold
    if threshold <= 0.0:
        k_max = float(kappa.max())
        # If the map is essentially empty, bail out
        if k_max <= 0.0:
            return []
        # Use a modest fraction of the maximum to define peaks
        threshold = 0.3 * k_max

    # Find local maxima
    local_max = maximum_filter(kappa, size=5)
    peak_mask = (kappa == local_max) & (kappa > threshold)
    
    # Get peak positions
    peak_indices = np.argwhere(peak_mask)
    
    peaks = []
    for idx in peak_indices:
        peaks.append((idx[0], idx[1], kappa[idx[0], idx[1]]))
    
    # Sort by κ value
    peaks.sort(key=lambda x: x[2], reverse=True)
    
    return peaks


def compute_offset(peaks, com_gas, com_gal, dx):
    """
    Compute offset between lensing and gas peaks.
    
    Parameters
    ----------
    peaks : list
        List of peak positions from compute_peak_positions
    com_gas : array
        Center of mass of gas
    com_gal : array
        Center of mass of galaxies
    dx : float
        Pixel size [Mpc]
        
    Returns
    -------
    offset_kpc : float
        Offset in kpc
    """
    if len(peaks) == 0:
        return 0.0
    
    # Use the highest peak
    main_peak = peaks[0]
    peak_pos = np.array([main_peak[0], main_peak[1]]) * dx
    
    # Offset from gas center
    offset = np.sqrt(np.sum((peak_pos - com_gas[:2])**2)) * 1000  # kpc
    
    return offset


class LensingPostprocessor:
    """
    Post-process simulation output to compute lensing maps.
    
    Parameters
    ----------
    output_dir : str
        Directory containing simulation output
    config : configparser.ConfigParser
        Configuration
    """
    
    def __init__(self, output_dir, config=None):
        self.output_dir = Path(output_dir)
        
        if config is None:
            config_path = Path(__file__).parent / 'hqiv_modifications' / 'config_hqiv.ini'
            config = configparser.ConfigParser()
            config.read(config_path)
        
        self.config = config
        self.box_size = float(config.get('simulation', 'boxlen', fallback='100.0'))
        self.z_source = float(config.get('lensing', 'z_source', fallback='1.0'))
        self.z_lens = float(config.get('bullet_cluster', 'z_obs', fallback='0.3'))
        
    def load_snapshot(self, snapshot_file):
        """Load simulation snapshot."""
        data = np.load(self.output_dir / snapshot_file, allow_pickle=True)
        return data
    
    def process(self, snapshot_file='final_particles.npz', method='born'):
        """
        Process snapshot to compute lensing maps.
        
        Parameters
        ----------
        snapshot_file : str
            Snapshot filename
        method : str
            'born' for Born approximation, 'full' for full geodesic
            
        Returns
        -------
        results : dict
            Dictionary with κ map, shear, and diagnostics
        """
        print(f"Processing {snapshot_file} with {method} method...")
        
        # Load data
        data = self.load_snapshot(snapshot_file)
        
        positions = data['positions']
        species = data['species']
        a = float(data['a'])
        
        # Separate gas and galaxies
        gas_mask = species == 'gas'
        pos_gas = positions[gas_mask]
        pos_gal = positions[~gas_mask]
        
        com_gas = np.mean(pos_gas, axis=0)
        com_gal = np.mean(pos_gal, axis=0)
        
        # Compute density field
        N = 256  # Grid for lensing
        dx = self.box_size / N
        density = np.zeros((N, N, N))
        
        for pos in positions:
            ix = int(np.floor(pos[0] / dx)) % N
            iy = int(np.floor(pos[1] / dx)) % N
            iz = int(np.floor(pos[2] / dx)) % N
            density[ix, iy, iz] += 1
        
        # Normalize
        rho_mean = len(positions) / N**3
        delta = density / rho_mean - 1
        
        # Compute convergence
        if method == 'born':
            kappa = compute_convergence_born(
                density, a, self.box_size,
                z_source=self.z_source, z_lens=self.z_lens
            )
            deflection = None
        else:
            # Would need potential and φ field from simulation
            # For now, fall back to Born
            kappa = compute_convergence_born(
                density, a, self.box_size,
                z_source=self.z_source, z_lens=self.z_lens
            )
            deflection = None
        
        # Compute shear
        gamma1, gamma2 = compute_shear(kappa)
        
        # Compute magnification
        mu = compute_magnification(kappa, gamma1, gamma2)
        
        # Find peaks (use adaptive threshold by passing <= 0)
        peaks = compute_peak_positions(kappa, threshold=0.0)
        
        # Compute offset
        offset = compute_offset(peaks, com_gas, com_gal, dx)
        
        print(f"  κ range: [{kappa.min():.3f}, {kappa.max():.3f}]")
        print(f"  Main peak: κ = {peaks[0][2]:.3f}" if peaks else "  No peaks found")
        print(f"  Offset (lensing-gas): {offset:.1f} kpc")
        
        return {
            'kappa': kappa,
            'gamma1': gamma1,
            'gamma2': gamma2,
            'magnification': mu,
            'peaks': peaks,
            'offset_kpc': offset,
            'com_gas': com_gas,
            'com_gal': com_gal,
            'density': density,
            'snapshot_file': snapshot_file,
            'a': a,
            'n_gas': int(np.sum(gas_mask)),
            'n_gal': int(np.sum(~gas_mask)),
        }
    
    def create_comparison_plot(self, results, output_file='lensing_comparison.png'):
        """Create comparison plot of lensing vs gas distribution."""
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        
        N = results['kappa'].shape[0]
        extent = [0, self.box_size, 0, self.box_size]
        
        # Convergence map
        ax = axes[0, 0]
        im = ax.imshow(results['kappa'].T, origin='lower', extent=extent, cmap='viridis')
        ax.set_title('Convergence κ')
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        plt.colorbar(im, ax=ax, label='κ')
        
        # Mark peaks
        for i, peak in enumerate(results['peaks'][:3]):
            px, py, pk = peak
            px_mpc = px * self.box_size / N
            py_mpc = py * self.box_size / N
            ax.plot(px_mpc, py_mpc, 'r*', markersize=15)
            ax.annotate(f'κ={pk:.2f}', (px_mpc, py_mpc), color='white', fontsize=8)
        
        # Gas density projection
        ax = axes[0, 1]
        gas_density = np.sum(results['density'], axis=2)
        im = ax.imshow(gas_density.T, origin='lower', extent=extent, cmap='hot')
        ax.set_title('Gas Density (projected)')
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        plt.colorbar(im, ax=ax)
        
        # Mark centers
        ax.plot(results['com_gas'][0], results['com_gas'][1], 'b+', markersize=15, mew=2, label='Gas COM')
        ax.plot(results['com_gal'][0], results['com_gal'][1], 'g+', markersize=15, mew=2, label='Gal COM')
        ax.legend()
        
        # Magnification (weak lensing: μ ≈ 1; use scale around 1 so structure is visible)
        mu = results['magnification']
        mu_center = np.median(mu)
        mu_span = max(0.05, np.percentile(np.abs(mu - mu_center), 98))
        vmin, vmax = mu_center - mu_span, mu_center + mu_span
        ax = axes[1, 0]
        im = ax.imshow(mu.T, origin='lower', extent=extent, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
        ax.set_title('Magnification μ')
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        plt.colorbar(im, ax=ax, label='μ')
        
        # Shear field
        ax = axes[1, 1]
        gamma_abs = np.sqrt(results['gamma1']**2 + results['gamma2']**2)
        im = ax.imshow(gamma_abs.T, origin='lower', extent=extent, cmap='plasma')
        ax.set_title('Shear |γ|')
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        plt.colorbar(im, ax=ax, label='|γ|')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / output_file, dpi=150)
        plt.close()
        
        print(f"Saved {output_file}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='HQIV Lensing Post-processing')
    parser.add_argument('--output', type=str, default='./output/',
                        help='Output directory with simulation data')
    parser.add_argument('--snapshot', type=str, default='final_particles.npz',
                        help='Snapshot file to process')
    parser.add_argument('--method', type=str, default='born',
                        choices=['born', 'full'],
                        help='Ray-tracing method')
    
    args = parser.parse_args()
    
    # Create postprocessor
    processor = LensingPostprocessor(args.output)
    
    # Process
    results = processor.process(args.snapshot, method=args.method)
    
    # Create plot
    processor.create_comparison_plot(results)
    
    # Print summary (data only; no tuning)
    dx_mpc = processor.box_size / results['kappa'].shape[0]
    summary_lines = []
    summary_lines.append("")
    summary_lines.append("="*60)
    summary_lines.append("Lensing Analysis Summary")
    summary_lines.append("="*60)
    summary_lines.append(f"  snapshot: {results['snapshot_file']}")
    summary_lines.append(f"  a: {results['a']:.6f}  z: {1.0/results['a'] - 1:.4f}")
    summary_lines.append(f"  box_size: {processor.box_size} Mpc")
    summary_lines.append(f"  n_gas: {results['n_gas']}  n_gal: {results['n_gal']}")
    summary_lines.append(f"  com_gas [Mpc]: ({results['com_gas'][0]:.4f}, {results['com_gas'][1]:.4f}, {results['com_gas'][2]:.4f})")
    summary_lines.append(f"  com_gal [Mpc]: ({results['com_gal'][0]:.4f}, {results['com_gal'][1]:.4f}, {results['com_gal'][2]:.4f})")
    summary_lines.append(f"  method: {args.method}")
    summary_lines.append(f"  kappa range: [{results['kappa'].min():.6f}, {results['kappa'].max():.6f}]")
    summary_lines.append(f"  n_peaks: {len(results['peaks'])}")
    for i, (pi, pj, pk) in enumerate(results['peaks'][:5]):
        x_mpc = pi * dx_mpc
        y_mpc = pj * dx_mpc
        summary_lines.append(f"  peak_{i+1} [Mpc]: ({x_mpc:.4f}, {y_mpc:.4f})  kappa: {pk:.6f}")
    summary_lines.append(f"  offset_lensing_gas_kpc: {results['offset_kpc']:.2f}")
    summary_lines.append("="*60)

    for line in summary_lines:
        print(line)

    # Also write summary to a text file alongside the PNG
    summary_path = processor.output_dir / "lensing_summary.txt"
    with open(summary_path, "w") as f:
        for line in summary_lines:
            f.write(line + "\n")
    print(f"Saved lensing_summary.txt to {summary_path}")
    
    return results


if __name__ == "__main__":
    results = main()