"""
HQIV Cosmology — Observables Module.

This module computes cosmological observables from the perturbation solutions:
- CMB angular power spectrum C_ℓ^TT
- Matter power spectrum P(k)
- Transfer functions
- σ₈

Paper Reference: Testing predictions against observations
- First acoustic peak position (standard: ℓ ~ 220-280)
- Growth suppression factor (~0.36× ΛCDM)
- Vorticity amplification

Author: HQIV Team
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import interp1d
from scipy.special import spherical_jn

# Import from other modules
from .background import c, Mpc_m, G0_SI, Gyr_s

# =============================================================================
# CMB Angular Power Spectrum
# =============================================================================

class CMBPowerSpectrum:
    """
    Compute CMB temperature anisotropy power spectrum C_ℓ^TT.
    
    Uses the line-of-sight integration method (Seljak & Zaldarriaga 1996)
    with simplified source functions for the scalar perturbations.
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology
    l_max : int
        Maximum multipole (default: 500)
    
    Notes
    -----
    The temperature anisotropy is:
        ΔT/T (n) = ∫ S(k, η) e^{ik·n(η₀-η)} dη
    
    where S(k, η) is the source function containing:
    - Sachs-Wolfe effect: Φ/3 at last scattering
    - Doppler term: v·n
    - Integrated Sachs-Wolfe: ∫ (Φ' + Ψ') dη
    - Acoustic oscillations: (δ_γ/4 + Φ)
    """
    
    def __init__(self, background, l_max=500):
        self.bg = background
        self.l_max = l_max
        
        # Recombination parameters
        self.z_rec = 1100  # Recombination redshift
        self.a_rec = 1.0 / (1.0 + self.z_rec)
        
        # Storage
        self.ell_arr = None
        self.Cl_TT = None
        self.Dl_TT = None
        
    def visibility_function(self, a):
        """
        Visibility function g(a) = e^{-τ} dτ/da.
        
        Peaks at recombination (a ~ 1/1100) with width determined by
        the thickness of the last scattering surface.
        
        Simplified model: Gaussian centered at recombination.
        """
        # Width of last scattering surface (Δz ~ 100)
        delta_a = 0.01 * self.a_rec
        
        # Gaussian visibility
        g = np.exp(-0.5 * ((a - self.a_rec) / delta_a)**2)
        
        # Normalize
        norm = np.sqrt(2 * np.pi) * delta_a
        return g / norm
    
    def source_function(self, k, a, delta, theta, Phi):
        """
        Compute the CMB source function S(k, η).
        
        Components:
        1. Sachs-Wolfe: Φ/3
        2. Acoustic: δ_γ/4 + Φ (approximated from matter δ)
        3. Doppler: v_b · n (from θ)
        4. ISW: Φ' (integrated along line of sight)
        
        Parameters
        ----------
        k : float
            Wavenumber [1/Mpc]
        a : float
            Scale factor
        delta : float
            Matter density contrast
        theta : float
            Velocity divergence
        Phi : float
            Newtonian potential
            
        Returns
        -------
        S : float
            Source function value
        """
        # Visibility function
        g = self.visibility_function(a)
        
        # Baryon-photon coupling at recombination
        # In the tight-coupling regime, δ_γ ≈ (4/3) δ_b
        
        # Sachs-Wolfe term
        SW = Phi / 3.0
        
        # Acoustic term (simplified)
        # The photon density contrast oscillates as cos(k r_s)
        # where r_s is the sound horizon
        eta = self.bg.get_eta(a) / Mpc_m  # [Mpc]
        r_s = self.bg.sound_horizon()  # [Mpc]
        
        # Acoustic oscillation
        acoustic = 0.5 * delta * np.cos(k * r_s * a / self.a_rec)
        
        # Doppler term (from velocity)
        # v ~ θ / k, and the Doppler contribution is v · n
        doppler = theta / (k + 0.01) * np.sin(k * r_s * a / self.a_rec)
        
        # ISW term (simplified: only late-time ISW)
        # Full treatment would integrate Φ' along line of sight
        ISW = 0.0  # Computed separately in line-of-sight integration
        
        # Total source
        S = g * (SW + acoustic + doppler) + ISW
        
        return S
    
    def compute_transfer_l(self, k, ell, eta_0, source_interp, eta_arr):
        """
        Compute the transfer function Δ_ℓ(k) for a single ℓ.
        
        Δ_ℓ(k) = ∫ S(k, η) j_ℓ(k(η₀ - η)) dη
        
        where j_ℓ is the spherical Bessel function.
        
        Parameters
        ----------
        k : float
            Wavenumber [1/Mpc]
        ell : int
            Multipole
        eta_0 : float
            Conformal time today [Mpc]
        source_interp : callable
            Interpolated source function S(η)
        eta_arr : array
            Conformal time array [Mpc]
            
        Returns
        -------
        Delta_l : float
            Transfer function value
        """
        # Integrate over conformal time
        integral = 0.0
        
        for i in range(1, len(eta_arr)):
            eta = eta_arr[i]
            eta_prev = eta_arr[i-1]
            deta = eta - eta_prev
            
            # Source function at this time
            S = source_interp(eta)
            
            # Spherical Bessel function argument
            x = k * (eta_0 - eta)
            
            if x > 0 and x < 1000:  # Avoid numerical issues
                # Use scipy's spherical Bessel function
                j_l = spherical_jn(ell, x, derivative=False)
            else:
                j_l = 0.0
            
            # Trapezoidal integration
            S_prev = source_interp(eta_prev)
            x_prev = k * (eta_0 - eta_prev)
            if x_prev > 0 and x_prev < 1000:
                j_l_prev = spherical_jn(ell, x_prev, derivative=False)
            else:
                j_l_prev = 0.0
            
            integral += 0.5 * (S * j_l + S_prev * j_l_prev) * deta
        
        return integral
    
    def compute_Cl(self, k_array, delta_k_arr, theta_k_arr, Phi_k_arr, 
                   a_start=1e-9, a_end=1.0, n_time=200):
        """
        Compute C_ℓ^TT from perturbation solutions.
        
        This is a simplified computation using the matter perturbations
        as a proxy for the full photon-baryon fluid evolution.
        
        Parameters
        ----------
        k_array : array
            Wavenumbers [1/Mpc]
        delta_k_arr : array
            Density contrast at recombination for each k
        theta_k_arr : array
            Velocity divergence at recombination for each k
        Phi_k_arr : array
            Potential at recombination for each k
        a_start : float
            Initial scale factor
        a_end : float
            Final scale factor
        n_time : int
            Number of time steps for line-of-sight integration
            
        Returns
        -------
        ell_arr : array
            Multipole array
        Cl_TT : array
            Temperature power spectrum C_ℓ
        Dl_TT : array
            D_ℓ = ℓ(ℓ+1)C_ℓ / 2π [μK²]
        """
        # Ensure background is computed
        if self.bg._eta_arr is None:
            self.bg.compute_background()
        
        # Conformal time array
        a_arr = np.logspace(np.log10(a_start), np.log10(a_end), n_time)
        eta_arr = np.array([self.bg.get_eta(a) / Mpc_m for a in a_arr])  # [Mpc]
        eta_0 = eta_arr[-1]  # Today
        
        # Primordial power spectrum (scale-invariant for simplicity)
        # P_Φ(k) = A_s * (k/k_pivot)^(n_s - 1)
        A_s = 2.1e-9  # Amplitude
        k_pivot = 0.05  # Mpc⁻¹
        n_s = 0.965  # Spectral index
        
        P_phi = A_s * (k_array / k_pivot)**(n_s - 1)
        
        # Multipole array
        self.ell_arr = np.arange(2, self.l_max + 1)
        self.Cl_TT = np.zeros(len(self.ell_arr))
        
        # Compute transfer functions for each ℓ
        print("Computing CMB transfer functions...")
        
        for i, ell in enumerate(self.ell_arr):
            if (i + 1) % 50 == 0:
                print(f"  ℓ = {ell} ({i+1}/{len(self.ell_arr)})")
            
            # Sum over k modes
            Cl_sum = 0.0
            
            for j, k in enumerate(k_array):
                # Create source function interpolator for this k
                # Simplified: use values at recombination
                delta_rec = delta_k_arr[j]
                theta_rec = theta_k_arr[j]
                Phi_rec = Phi_k_arr[j]
                
                # Source function (simplified)
                def source_func(eta):
                    # Find corresponding a
                    a = np.interp(eta, eta_arr, a_arr)
                    g = self.visibility_function(a)
                    # Simplified source
                    return g * (Phi_rec / 3.0 + 0.5 * delta_rec)
                
                # Transfer function
                Delta_l = self.compute_transfer_l(
                    k, ell, eta_0, 
                    interp1d(eta_arr, [source_func(e) for e in eta_arr], 
                            kind='linear', fill_value=0.0, bounds_error=False),
                    eta_arr
                )
                
                # Add to power spectrum
                # C_ℓ = (2/π) ∫ k² P_Φ(k) |Δ_ℓ(k)|² dk
                dk = k_array[1] - k_array[0] if j > 0 else k_array[1] - k_array[0]
                Cl_sum += k**2 * P_phi[j] * Delta_l**2 * dk
            
            self.Cl_TT[i] = 2.0 / np.pi * Cl_sum
        
        # Convert to D_ℓ = ℓ(ℓ+1) C_ℓ / 2π in μK²
        # CMB temperature
        T_CMB_uK = 2.7255e6  # μK
        
        self.Dl_TT = self.ell_arr * (self.ell_arr + 1) * self.Cl_TT * T_CMB_uK**2 / (2 * np.pi)
        
        return self.ell_arr, self.Cl_TT, self.Dl_TT
    
    def find_acoustic_peak(self):
        """
        Find the position of the first acoustic peak.
        
        Standard ΛCDM: ℓ_peak ~ 220-280
        
        Returns
        -------
        ell_peak : float
            Position of first acoustic peak
        """
        if self.Dl_TT is None:
            raise RuntimeError("Run compute_Cl() first")
        
        # Find maximum in D_ℓ
        # Look in range ℓ = 100-400
        mask = (self.ell_arr >= 100) & (self.ell_arr <= 400)
        if not np.any(mask):
            return None
        
        ell_subset = self.ell_arr[mask]
        Dl_subset = self.Dl_TT[mask]
        
        # Find peak
        peak_idx = np.argmax(Dl_subset)
        ell_peak = ell_subset[peak_idx]
        
        return ell_peak


def compute_matter_power_spectrum(background, k_array, delta_k_arr, a=1.0):
    """
    Compute the matter power spectrum P(k).
    
    P(k) = ⟨|δ(k)|²⟩ = A_s × T(k)² × D(a)²
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology
    k_array : array
        Wavenumbers [1/Mpc]
    delta_k_arr : array
        Density contrast at scale factor a for each k
    a : float
        Scale factor (default: a = 1, today)
        
    Returns
    -------
    P_k : array
        Matter power spectrum [Mpc³]
    """
    # Primordial power spectrum
    A_s = 2.1e-9
    k_pivot = 0.05  # Mpc⁻¹
    n_s = 0.965
    
    P_primordial = A_s * (k_array / k_pivot)**(n_s - 1)
    
    # Transfer function from δ
    T_k = delta_k_arr / delta_k_arr[0] if delta_k_arr[0] != 0 else np.ones_like(delta_k_arr)
    
    # Normalize
    T_k = T_k / np.max(np.abs(T_k))
    
    # Matter power spectrum
    P_k = P_primordial * T_k**2
    
    # Scale factor normalization
    # P(k, a) = P(k) × D(a)²
    # For a = 1, D = 1
    
    return P_k


def compute_sigma8(background, k_array, P_k):
    """
    Compute σ₈, the RMS density fluctuation in 8 Mpc/h spheres.
    
    σ₈² = (1/2π²) ∫ k³ P(k) W²(kR₈) dk
    
    where W(kR) = 3 × (sin(kR) - kR cos(kR)) / (kR)³
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology
    k_array : array
        Wavenumbers [1/Mpc]
    P_k : array
        Matter power spectrum [Mpc³]
        
    Returns
    -------
    sigma8 : float
        RMS density fluctuation
    """
    # R₈ = 8 Mpc/h
    h = background.H0_km / 100.0
    R8 = 8.0 / h  # Mpc
    
    # Top-hat window function in Fourier space
    def W_kR(kR):
        return 3.0 * (np.sin(kR) - kR * np.cos(kR)) / (kR**3 + 1e-10)
    
    # Integrate
    integral = 0.0
    for i in range(1, len(k_array)):
        k = k_array[i]
        k_prev = k_array[i-1]
        dk = k - k_prev
        
        kR = k * R8
        kR_prev = k_prev * R8
        
        W2 = W_kR(kR)**2
        W2_prev = W_kR(kR_prev)**2
        
        integrand = k**3 * P_k[i] * W2
        integrand_prev = k_prev**3 * P_k[i-1] * W2_prev
        
        integral += 0.5 * (integrand + integrand_prev) * dk
    
    sigma8_sq = integral / (2 * np.pi**2)
    
    return np.sqrt(sigma8_sq)


# =============================================================================
# Simplified CMB Computation (for quick testing)
# =============================================================================

def simplified_CMB_spectrum(background, l_max=500, n_k=50):
    """
    Compute a simplified CMB spectrum for quick testing.
    
    Uses analytic approximations for the acoustic oscillations
    without full line-of-sight integration.
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology
    l_max : int
        Maximum multipole
    n_k : int
        Number of k modes
        
    Returns
    -------
    ell_arr : array
        Multipole array
    Dl_TT : array
        D_ℓ spectrum [μK²]
    ell_peak : float
        First acoustic peak position
    """
    # Ensure background is computed
    if background._eta_arr is None:
        background.compute_background()
    
    # Key scales
    eta_0 = background.conformal_time_today()  # Mpc
    r_s = background.sound_horizon()  # Mpc
    
    # Angular diameter distance to recombination
    a_rec = 1.0 / 1100.0
    z_rec = 1100.0
    
    # d_A to recombination
    # In flat universe: d_A = η_0 - η_rec
    eta_rec = background.get_eta(a_rec) / Mpc_m
    d_A = (eta_0 - eta_rec)  # Mpc
    
    # Peak position: ℓ_peak ≈ π d_A / r_s
    ell_peak_theory = np.pi * d_A / r_s
    
    print(f"Sound horizon: {r_s:.2f} Mpc")
    print(f"Angular diameter distance: {d_A:.2f} Mpc")
    print(f"Predicted peak position: ℓ ≈ {ell_peak_theory:.0f}")
    
    # Multipole array
    ell_arr = np.arange(2, l_max + 1, dtype=float)
    
    # Simplified spectrum model
    # D_ℓ = A × (1 + 7 × 10⁻⁶ × ℓ) × exp(-(ℓ/1500)²) × [SW + acoustic]
    
    # Sachs-Wolfe plateau (low ℓ)
    SW = 1000.0 * np.ones_like(ell_arr)  # μK²
    
    # Acoustic oscillations
    # Position of peaks: ℓ_n = n × π × d_A / r_s
    # First peak: n = 1
    # Second peak: n = 2, etc.
    
    acoustic = np.zeros_like(ell_arr)
    for n in range(1, 5):  # First 4 peaks
        ell_n = n * ell_peak_theory
        amplitude = 6000.0 / n**1.5  # Decreasing amplitude
        width = 50.0 + 20.0 * n
        acoustic += amplitude * np.exp(-0.5 * ((ell_arr - ell_n) / width)**2)
    
    # Damping at high ℓ (diffusion damping)
    damping = np.exp(-(ell_arr / 1500.0)**2)
    
    # Total spectrum
    Dl_TT = (SW + acoustic) * damping
    
    # Find actual peak
    mask = (ell_arr >= 100) & (ell_arr <= 400)
    peak_idx = np.argmax(Dl_TT[mask])
    ell_peak = ell_arr[mask][peak_idx]
    
    return ell_arr, Dl_TT, ell_peak


def compare_CMB_peaks(hqiv_bg, lcdm_bg, l_max=500):
    """
    Compare CMB peak positions between HQIV and ΛCDM.
    
    Parameters
    ----------
    hqiv_bg : CosmologicalBackground
        HQIV background
    lcdm_bg : CosmologicalBackground
        ΛCDM background
    l_max : int
        Maximum multipole
        
    Returns
    -------
    comparison : dict
        Dictionary with peak positions and ratios
    """
    print("\nComparing CMB peak positions...")
    
    # Compute spectra
    ell_h, Dl_h, peak_h = simplified_CMB_spectrum(hqiv_bg, l_max)
    ell_l, Dl_l, peak_l = simplified_CMB_spectrum(lcdm_bg, l_max)
    
    # Results
    comparison = {
        'ell_HQIV': peak_h,
        'ell_LCDM': peak_l,
        'ratio': peak_h / peak_l,
        'shift': peak_h - peak_l,
        'ell_arr': ell_h,
        'Dl_HQIV': Dl_h,
        'Dl_LCDM': Dl_l
    }
    
    print(f"\nHQIV first peak: ℓ = {peak_h:.0f}")
    print(f"ΛCDM first peak: ℓ = {peak_l:.0f}")
    print(f"Shift: Δℓ = {peak_h - peak_l:.0f}")
    
    return comparison


# =============================================================================
# Testing and visualization
# =============================================================================

def test_CMB_spectrum():
    """Test CMB spectrum computation and create plots."""
    import matplotlib.pyplot as plt
    
    print("Testing CMB spectrum computation...")
    
    # Create background cosmologies
    hqiv_bg = CosmologicalBackground(H0=73.2, Om_m=0.031, hqiv_on=True,
                                       beta=1.02, Om_h=1.00, n_h=1.04)
    lcdm_bg = CosmologicalBackground(H0=67.4, Om_m=0.315, hqiv_on=False)
    
    hqiv_bg.compute_background()
    lcdm_bg.compute_background()
    
    # Compare peaks
    comparison = compare_CMB_peaks(hqiv_bg, lcdm_bg, l_max=500)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    ell_arr = comparison['ell_arr']
    
    # CMB spectrum comparison
    ax = axes[0, 0]
    ax.plot(ell_arr, comparison['Dl_HQIV'], 'b-', label='HQIV', linewidth=2)
    ax.plot(ell_arr, comparison['Dl_LCDM'], 'r--', label='ΛCDM', linewidth=2)
    ax.axvline(comparison['ell_HQIV'], color='b', linestyle=':', alpha=0.5)
    ax.axvline(comparison['ell_LCDM'], color='r', linestyle=':', alpha=0.5)
    ax.set_xlabel('Multipole ℓ')
    ax.set_ylabel('D_ℓ^TT [μK²]')
    ax.set_title('CMB Temperature Power Spectrum')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim([2, 500])
    
    # Ratio
    ax = axes[0, 1]
    ratio = comparison['Dl_HQIV'] / (comparison['Dl_LCDM'] + 1e-10)
    ax.plot(ell_arr, ratio, 'b-', linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Multipole ℓ')
    ax.set_ylabel('D_ℓ^HQIV / D_ℓ^ΛCDM')
    ax.set_title('Spectrum Ratio')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([2, 500])
    
    # Sound horizon comparison
    ax = axes[1, 0]
    r_s_h = hqiv_bg.sound_horizon()
    r_s_l = lcdm_bg.sound_horizon()
    d_A_h = hqiv_bg.conformal_time_today() - hqiv_bg.get_eta(1/1100) / Mpc_m
    d_A_l = lcdm_bg.conformal_time_today() - lcdm_bg.get_eta(1/1100) / Mpc_m
    
    labels = ['Sound horizon\nr_s [Mpc]', 'Ang. diam. dist.\nd_A [Mpc]', 'Peak pos.\nℓ = πd_A/r_s']
    hqiv_vals = [r_s_h, d_A_h, np.pi * d_A_h / r_s_h]
    lcdm_vals = [r_s_l, d_A_l, np.pi * d_A_l / r_s_l]
    
    x = np.arange(len(labels))
    width = 0.35
    
    ax.bar(x - width/2, hqiv_vals, width, label='HQIV', color='blue', alpha=0.7)
    ax.bar(x + width/2, lcdm_vals, width, label='ΛCDM', color='red', alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Value')
    ax.set_title('Key CMB Scales')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Peak position summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_text = f"""
    CMB First Acoustic Peak Position
    {'='*40}
    
    HQIV:
      Sound horizon: {r_s_h:.1f} Mpc
      d_A to recombination: {d_A_h:.1f} Mpc
      Predicted peak: ℓ ≈ {np.pi * d_A_h / r_s_h:.0f}
    
    ΛCDM:
      Sound horizon: {r_s_l:.1f} Mpc
      d_A to recombination: {d_A_l:.1f} Mpc
      Predicted peak: ℓ ≈ {np.pi * d_A_l / r_s_l:.0f}
    
    Peak shift: Δℓ = {comparison['shift']:.0f}
    Ratio: {comparison['ratio']:.3f}
    
    Standard ΛCDM peak: ℓ ~ 220-280
    """
    
    ax.text(0.1, 0.5, summary_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='center', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('hqiv_solver/cmb_spectrum.png', dpi=150)
    plt.close()
    print("\nSaved cmb_spectrum.png")
    
    return comparison


if __name__ == "__main__":
    test_CMB_spectrum()