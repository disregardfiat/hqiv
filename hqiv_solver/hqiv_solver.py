"""
HQIV Cosmology — Main Solver Class.

This module provides the main HQIVPerturbations class that orchestrates
the full perturbation evolution including scalar, vector, and tensor modes.

Paper Reference: Full framework implementation
- Background: Section 5
- Scalar perturbations: Section 7, Eqs. 9-11
- Vector perturbations: Section 7, Eq. 12
- Observables: Section 8

Author: HQIV Team
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import time

# Import modules
from .background import CosmologicalBackground, c, Mpc_m
from .scalars import ScalarPerturbations, compute_transfer_function, inertia_reduction_factor
from .vectors import VorticityPerturbations, VorticityScalarCoupling, compute_vorticity_growth_exponent
from .observables import CMBPowerSpectrum, simplified_CMB_spectrum, compute_sigma8


class HQIVPerturbations:
    """
    Master class for HQIV perturbation evolution.
    
    This class provides a unified interface to:
    1. Set up background cosmology (ΛCDM or HQIV)
    2. Evolve scalar perturbations with optional HQIV modifications
    3. Evolve vector (vorticity) perturbations with horizon coupling
    4. Compute observables: C_ℓ, P(k), σ₈, growth factor
    
    Parameters
    ----------
    H0 : float
        Hubble constant today [km/s/Mpc]
    Om_m : float
        Matter density parameter today
    Om_r : float, optional
        Radiation density parameter. If None, computed from CMB temperature.
    hqiv_on : bool
        If True, use HQIV modifications. If False, use standard ΛCDM.
    beta : float
        HQIV horizon-smoothing parameter (paper eq. 1)
    Om_h : float
        HQIV effective horizon density parameter
    n_h : float
        HQIV horizon dilution exponent
    alpha_G : float
        Exponent for varying G(a)
    inertia_form : str
        Form of inertia reduction: 'sqrt' or 'brodie'
    
    Example
    -------
    >>> solver = HQIVPerturbations(H0=73.2, Om_m=0.031, hqiv_on=True)
    >>> solver.run(k_array=np.logspace(-3, 0, 50))
    >>> solver.plot_results()
    """
    
    def __init__(self, H0=73.2, Om_m=0.031, Om_r=None, hqiv_on=True,
                 beta=1.02, Om_h=1.00, n_h=1.04, alpha_G=0.6, 
                 inertia_form='sqrt'):
        
        # Store parameters
        self.params = {
            'H0': H0,
            'Om_m': Om_m,
            'Om_r': Om_r,
            'hqiv_on': hqiv_on,
            'beta': beta,
            'Om_h': Om_h,
            'n_h': n_h,
            'alpha_G': alpha_G,
            'inertia_form': inertia_form
        }
        
        # Create background
        self.background = CosmologicalBackground(
            H0=H0, Om_m=Om_m, Om_r=Om_r,
            hqiv_on=hqiv_on, beta=beta, Om_h=Om_h, 
            n_h=n_h, alpha_G=alpha_G
        )
        
        # Storage for results
        self.k_array = None
        self.scalar_results = {}
        self.vector_results = {}
        self.observables = {}
        self.comparison = {}
        
    def run(self, k_array=None, a_start=1e-9, a_end=1.0, n_pts=500,
            compute_vectors=True, compute_observables=True, verbose=True):
        """
        Run the full perturbation evolution for all k modes.
        
        Parameters
        ----------
        k_array : array, optional
            Wavenumbers [1/Mpc]. If None, uses default range.
        a_start : float
            Initial scale factor
        a_end : float
            Final scale factor
        n_pts : int
            Number of time steps
        compute_vectors : bool
            Whether to compute vorticity evolution
        compute_observables : bool
            Whether to compute CMB and other observables
        verbose : bool
            Print progress information
            
        Returns
        -------
        results : dict
            Dictionary containing all computed quantities
        """
        start_time = time.time()
        
        # Default k array
        if k_array is None:
            k_array = np.logspace(-3, 0, 50)  # 0.001 to 1 Mpc⁻¹
        
        self.k_array = k_array
        
        if verbose:
            print("=" * 60)
            print("HQIV Perturbation Solver")
            print("=" * 60)
            print(f"Mode: {'HQIV' if self.params['hqiv_on'] else 'ΛCDM'}")
            print(f"k range: {k_array[0]:.4f} - {k_array[-1]:.4f} Mpc⁻¹")
            print(f"a range: {a_start:.2e} - {a_end:.2e}")
            print("-" * 60)
        
        # Step 1: Compute background
        if verbose:
            print("\n[1/4] Computing background cosmology...")
        
        self.background.compute_background(a_start, a_end, n_pts)
        
        if verbose:
            self.background.summary()
        
        # Step 2: Compute scalar perturbations for each k
        if verbose:
            print("\n[2/4] Computing scalar perturbations...")
        
        delta_k = np.zeros(len(k_array))
        theta_k = np.zeros(len(k_array))
        Phi_k = np.zeros(len(k_array))
        
        # Store evolution for a few representative k values
        k_representative = [0.01, 0.1, 0.5]  # Mpc⁻¹
        
        for i, k in enumerate(k_array):
            pert = ScalarPerturbations(
                self.background, k, 
                hqiv_on=self.params['hqiv_on'],
                inertia_form=self.params['inertia_form']
            )
            a_arr, delta_arr, theta_arr, Phi_arr = pert.solve(a_start, a_end, n_pts)
            
            delta_k[i] = delta_arr[-1]
            theta_k[i] = theta_arr[-1]
            Phi_k[i] = Phi_arr[-1]
            
            # Store full evolution for representative k
            if k in k_representative or i == len(k_array) // 2:
                self.scalar_results[f'k_{k:.3f}'] = {
                    'a': a_arr,
                    'delta': delta_arr,
                    'theta': theta_arr,
                    'Phi': Phi_arr
                }
            
            if verbose and (i + 1) % 10 == 0:
                print(f"  Progress: {i+1}/{len(k_array)} k modes")
        
        self.scalar_results['k_array'] = k_array
        self.scalar_results['delta_final'] = delta_k
        self.scalar_results['theta_final'] = theta_k
        self.scalar_results['Phi_final'] = Phi_k
        
        # Compute growth factor
        k_test = 0.1  # Mpc⁻¹
        pert_test = ScalarPerturbations(
            self.background, k_test,
            hqiv_on=self.params['hqiv_on']
        )
        a_arr, delta_arr, _, _ = pert_test.solve(a_start, a_end, n_pts)
        D_arr = delta_arr / delta_arr[-1]
        
        self.scalar_results['growth_factor'] = {
            'a': a_arr,
            'D': D_arr
        }
        
        # Step 3: Compute vector perturbations
        if compute_vectors:
            if verbose:
                print("\n[3/4] Computing vorticity evolution...")
            
            vort = VorticityPerturbations(
                self.background, k_test,
                hqiv_on=self.params['hqiv_on']
            )
            a_vort, omega_arr = vort.solve(a_start, a_end, n_pts)
            omega_mag = vort.vorticity_magnitude()
            
            self.vector_results = {
                'a': a_vort,
                'omega': omega_arr,
                'omega_mag': omega_mag,
                'amplification': omega_mag / omega_mag[0]
            }
            
            # Compute vorticity growth exponent
            beta_eff, _, _ = compute_vorticity_growth_exponent(
                self.background, k_test, a_start, a_end
            )
            self.vector_results['beta_eff'] = beta_eff
        else:
            self.vector_results = None
        
        # Step 4: Compute observables
        if compute_observables:
            if verbose:
                print("\n[4/4] Computing observables...")
            
            # CMB spectrum (simplified)
            ell_arr, Dl_TT, ell_peak = simplified_CMB_spectrum(
                self.background, l_max=500
            )
            
            self.observables = {
                'ell': ell_arr,
                'Dl_TT': Dl_TT,
                'ell_peak': ell_peak,
                'sound_horizon': self.background.sound_horizon(),
                'conformal_time': self.background.conformal_time_today()
            }
            
            # σ₈ (simplified)
            # Would need proper P(k) calculation
        else:
            self.observables = None
        
        elapsed = time.time() - start_time
        
        if verbose:
            print("\n" + "=" * 60)
            print(f"Computation complete in {elapsed:.2f} seconds")
            print("=" * 60)
        
        return {
            'background': self.background,
            'scalars': self.scalar_results,
            'vectors': self.vector_results,
            'observables': self.observables
        }
    
    def compare_with_lcdm(self, k_array=None, a_start=1e-9, a_end=1.0, 
                          n_pts=500, verbose=True):
        """
        Compare HQIV results with standard ΛCDM.
        
        Creates a ΛCDM cosmology with the same H0 and Om_m for comparison.
        
        Parameters
        ----------
        k_array : array, optional
            Wavenumbers [1/Mpc]
        a_start : float
            Initial scale factor
        a_end : float
            Final scale factor
        n_pts : int
            Number of time steps
        verbose : bool
            Print progress
            
        Returns
        -------
        comparison : dict
            Dictionary with comparison results
        """
        if verbose:
            print("\n" + "=" * 60)
            print("Comparing HQIV with ΛCDM")
            print("=" * 60)
        
        # Create ΛCDM comparison
        # Use standard ΛCDM parameters for fair comparison
        lcdm = HQIVPerturbations(
            H0=67.4,  # Planck value
            Om_m=0.315,  # Planck value
            hqiv_on=False
        )
        
        # Run ΛCDM
        if verbose:
            print("\nRunning ΛCDM...")
        lcdm_results = lcdm.run(
            k_array=k_array, a_start=a_start, a_end=a_end,
            n_pts=n_pts, compute_vectors=True, compute_observables=True,
            verbose=False
        )
        
        # Run HQIV if not already done
        if self.scalar_results is None or len(self.scalar_results) == 0:
            if verbose:
                print("\nRunning HQIV...")
            hqiv_results = self.run(
                k_array=k_array, a_start=a_start, a_end=a_end,
                n_pts=n_pts, compute_vectors=True, compute_observables=True,
                verbose=False
            )
        else:
            hqiv_results = {
                'background': self.background,
                'scalars': self.scalar_results,
                'vectors': self.vector_results,
                'observables': self.observables
            }
        
        # Compare growth factors
        D_hqiv = hqiv_results['scalars']['growth_factor']['D']
        a_hqiv = hqiv_results['scalars']['growth_factor']['a']
        D_lcdm = lcdm_results['scalars']['growth_factor']['D']
        a_lcdm = lcdm_results['scalars']['growth_factor']['a']
        
        # Interpolate to common grid
        a_common = a_hqiv
        D_lcdm_interp = np.interp(a_common, a_lcdm, D_lcdm)
        growth_suppression = D_hqiv / D_lcdm_interp
        
        # Compare vorticity
        if hqiv_results['vectors'] and lcdm_results['vectors']:
            vort_ratio = (hqiv_results['vectors']['amplification'][-1] / 
                         lcdm_results['vectors']['amplification'][-1])
        else:
            vort_ratio = None
        
        # Compare CMB peaks
        if hqiv_results['observables'] and lcdm_results['observables']:
            peak_shift = (hqiv_results['observables']['ell_peak'] - 
                         lcdm_results['observables']['ell_peak'])
        else:
            peak_shift = None
        
        self.comparison = {
            'a': a_common,
            'D_HQIV': D_hqiv,
            'D_LCDM': D_lcdm_interp,
            'growth_suppression': growth_suppression,
            'vorticity_ratio': vort_ratio,
            'peak_shift': peak_shift,
            'hqiv_results': hqiv_results,
            'lcdm_results': lcdm_results
        }
        
        if verbose:
            print("\n" + "-" * 60)
            print("Comparison Summary:")
            print("-" * 60)
            print(f"Growth suppression at a=0.5: {np.interp(0.5, a_common, growth_suppression):.4f}")
            print(f"  (Paper prediction: ~0.36)")
            if vort_ratio:
                print(f"Vorticity amplification ratio: {vort_ratio:.2e}")
            if peak_shift:
                print(f"CMB peak shift: Δℓ = {peak_shift:.0f}")
            print("-" * 60)
        
        return self.comparison
    
    def plot_results(self, save_path='hqiv_solver/results.png'):
        """
        Create summary plots of all results.
        
        Parameters
        ----------
        save_path : str
            Path to save the figure
        """
        import matplotlib.pyplot as plt
        
        if self.scalar_results is None or len(self.scalar_results) == 0:
            print("No results to plot. Run the solver first.")
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Plot 1: Background H(a)
        ax = axes[0, 0]
        a_bg = self.background._a_arr
        H_bg = self.background.H_km(a_bg)
        ax.loglog(a_bg, H_bg, 'b-', linewidth=2)
        ax.set_xlabel('Scale factor a')
        ax.set_ylabel('H(a) [km/s/Mpc]')
        ax.set_title('Hubble Parameter')
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Growth factor
        ax = axes[0, 1]
        if 'growth_factor' in self.scalar_results:
            a_D = self.scalar_results['growth_factor']['a']
            D = self.scalar_results['growth_factor']['D']
            ax.semilogx(a_D, D, 'b-', linewidth=2)
            ax.set_xlabel('Scale factor a')
            ax.set_ylabel('D(a)')
            ax.set_title('Linear Growth Factor')
            ax.grid(True, alpha=0.3)
        
        # Plot 3: Transfer function
        ax = axes[0, 2]
        if 'k_array' in self.scalar_results:
            k = self.scalar_results['k_array']
            delta = self.scalar_results['delta_final']
            T = delta / np.max(np.abs(delta))
            ax.semilogx(k, T, 'b-', linewidth=2)
            ax.set_xlabel('k [Mpc⁻¹]')
            ax.set_ylabel('T(k)')
            ax.set_title('Transfer Function')
            ax.grid(True, alpha=0.3)
        
        # Plot 4: Vorticity evolution
        ax = axes[1, 0]
        if self.vector_results:
            a_v = self.vector_results['a']
            A_omega = self.vector_results['amplification']
            ax.loglog(a_v, A_omega, 'b-', linewidth=2, label='HQIV')
            # Standard decay for comparison
            ax.loglog(a_v, (a_v / a_v[0])**(-2), 'r--', linewidth=1.5, 
                     label='Standard: a⁻²', alpha=0.7)
            ax.set_xlabel('Scale factor a')
            ax.set_ylabel('|ω| / |ω₀|')
            ax.set_title('Vorticity Evolution')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Plot 5: CMB spectrum
        ax = axes[1, 1]
        if self.observables:
            ell = self.observables['ell']
            Dl = self.observables['Dl_TT']
            ax.plot(ell, Dl, 'b-', linewidth=2)
            ax.axvline(self.observables['ell_peak'], color='r', 
                      linestyle=':', alpha=0.7, 
                      label=f"Peak: ℓ = {self.observables['ell_peak']:.0f}")
            ax.set_xlabel('Multipole ℓ')
            ax.set_ylabel('D_ℓ^TT [μK²]')
            ax.set_title('CMB Temperature Spectrum')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_xlim([2, 500])
        
        # Plot 6: Comparison with ΛCDM (if available)
        ax = axes[1, 2]
        if self.comparison:
            a_c = self.comparison['a']
            supp = self.comparison['growth_suppression']
            ax.semilogx(a_c, supp, 'b-', linewidth=2)
            ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
            ax.axhline(0.36, color='g', linestyle=':', alpha=0.7, 
                      label='Paper prediction (0.36)')
            ax.set_xlabel('Scale factor a')
            ax.set_ylabel('D_HQIV / D_ΛCDM')
            ax.set_title('Growth Suppression vs ΛCDM')
            ax.legend()
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'Run compare_with_lcdm()\nfor comparison',
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=12, color='gray')
            ax.axis('off')
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=150)
        plt.close()
        print(f"\nSaved results to {save_path}")
    
    def summary(self):
        """Print a summary of results."""
        print("\n" + "=" * 60)
        print("HQIV Perturbation Solver - Summary")
        print("=" * 60)
        
        print("\nParameters:")
        for key, val in self.params.items():
            print(f"  {key}: {val}")
        
        if self.background._a_arr is not None:
            print("\nBackground:")
            print(f"  Universe age: {self.background.age_today():.2f} Gyr")
            print(f"  Sound horizon: {self.background.sound_horizon():.2f} Mpc")
        
        if self.scalar_results and 'growth_factor' in self.scalar_results:
            print("\nScalar Perturbations:")
            D = self.scalar_results['growth_factor']['D']
            a = self.scalar_results['growth_factor']['a']
            print(f"  Growth factor at a=0.5: {np.interp(0.5, a, D):.4f}")
        
        if self.vector_results:
            print("\nVector Perturbations:")
            print(f"  Vorticity growth exponent: {self.vector_results['beta_eff']:.4f}")
            print(f"  Final amplification: {self.vector_results['amplification'][-1]:.2e}")
        
        if self.observables:
            print("\nObservables:")
            print(f"  CMB first peak: ℓ = {self.observables['ell_peak']:.0f}")
        
        if self.comparison:
            print("\nComparison with ΛCDM:")
            a_c = self.comparison['a']
            supp = self.comparison['growth_suppression']
            print(f"  Growth suppression at a=0.5: {np.interp(0.5, a_c, supp):.4f}")
            if self.comparison['peak_shift']:
                print(f"  CMB peak shift: Δℓ = {self.comparison['peak_shift']:.0f}")
        
        print("=" * 60)


# =============================================================================
# Convenience functions
# =============================================================================

def quick_test(hqiv_on=True, compare=True):
    """
    Run a quick test of the solver.
    
    Parameters
    ----------
    hqiv_on : bool
        Test HQIV or ΛCDM
    compare : bool
        Compare with the other model
    """
    print("Running quick test...")
    
    # Create solver
    solver = HQIVPerturbations(
        H0=73.2 if hqiv_on else 67.4,
        Om_m=0.031 if hqiv_on else 0.315,
        hqiv_on=hqiv_on
    )
    
    # Run with default parameters
    k_array = np.logspace(-2, 0, 20)  # Fewer modes for speed
    solver.run(k_array=k_array, n_pts=200, verbose=True)
    
    # Compare with ΛCDM
    if compare:
        solver.compare_with_lcdm(k_array=k_array, n_pts=200, verbose=True)
    
    # Plot results
    solver.plot_results()
    
    # Print summary
    solver.summary()
    
    return solver


if __name__ == "__main__":
    quick_test(hqiv_on=True, compare=True)