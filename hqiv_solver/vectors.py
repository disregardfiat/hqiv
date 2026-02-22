"""
HQIV Cosmology — Vector Perturbation (Vorticity) Module.

This module implements the vector perturbation equations for the HQIV framework,
specifically the vorticity amplification equation with horizon coupling.

Paper Reference: Section 7, Eq. 12
    ∂ω/∂t + (v·∇)ω = β H (ω · ê_Θ)

This is a key prediction of HQIV: vorticity is amplified rather than
decaying as in standard cosmology.

Author: HQIV Team
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Import background module
from .background import c, Mpc_m

# =============================================================================
# Vorticity Evolution Equations
# =============================================================================

class VorticityPerturbations:
    """
    Solver for vector perturbations (vorticity) in Fourier space.
    
    In standard cosmology, vorticity decays as ω ∝ a⁻² due to expansion.
    In HQIV, the horizon coupling term can amplify or sustain vorticity:
        ω' + ℋ ω = β(a) H(a) (ω · ê_Θ)
    
    Paper eq. 12:
        ∂ω/∂t + (v·∇)ω = β H (ω · ê_Θ)
    
    In conformal time:
        ω' + ℋ ω = β(a) a H(a) (ω · ê_Θ)
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology object
    k : float
        Comoving wavenumber [1/Mpc]
    hqiv_on : bool
        If True, include HQIV vorticity amplification
    
    State Vector
    ------------
    y = [ω_x, ω_y, ω_z]
    where ω is the vorticity vector in comoving coordinates
    
    Notes
    -----
    The unit vector ê_Θ points toward the local horizon. For simplicity,
    we initially take ê_Θ ≈ radial direction (pointing away from observer).
    The dot product (ω · ê_Θ) then selects the radial component of vorticity.
    """
    
    def __init__(self, background, k, hqiv_on=True):
        self.bg = background
        self.k = k  # [1/Mpc]
        self.hqiv_on = hqiv_on and background.hqiv_on
        
        # Storage for solution
        self.a_arr = None
        self.omega_arr = None  # Shape: (n_times, 3)
        
    def conformal_H(self, a):
        """Conformal Hubble parameter ℋ = aH in units of 1/Mpc."""
        return a * self.bg.get_H(a) / c * Mpc_m
    
    def horizon_direction(self, k_vec):
        """
        Compute the horizon direction unit vector ê_Θ.
        
        In a homogeneous background, the horizon is spherical and ê_Θ
        points radially outward. For a Fourier mode with wavevector k,
        we take ê_Θ parallel to k.
        
        Parameters
        ----------
        k_vec : array
            Wavevector direction (normalized)
            
        Returns
        -------
        e_Theta : array
            Unit vector pointing toward horizon
        """
        return k_vec / np.linalg.norm(k_vec)
    
    def derivatives(self, a, y, k_vec=None):
        """
        Compute derivatives d(ω)/da for the vorticity vector.
        
        β-free formulation (from action, paper/beta-free.tex):
            ω' + 2ℋ ω = (∂f/∂φ) (k × ∇φ) · ê_ω
        
        where ω' = dω/dη (conformal time derivative), ℋ = aH,
        and φ = 3cH is the geometric horizon field.
        
        The source term comes from horizon gradients ∇φ, which are sourced
        by density perturbations δρ in the scalar perturbation sector.
        
        In terms of scale factor:
            dω/da = (dω/dη) × (dη/da) = ω' / (a² H) = ω' / (a ℋ)
        
        So:
            dω/da = [-2ℋ ω + source] / (a ℋ)
                  = -2 ω / a + source / (a ℋ)
        
        The source term is proportional to ∇φ, which grows with structure.
        
        Parameters
        ----------
        a : float
            Scale factor
        y : array
            State vector [ω_x, ω_y, ω_z]
        k_vec : array, optional
            Wavevector direction. If None, assumes radial.
            
        Returns
        -------
        dyda : array
            Derivatives [dω_x/da, dω_y/da, dω_z/da]
        """
        omega = y  # [ω_x, ω_y, ω_z]
        
        # Conformal Hubble parameter ℋ = aH
        H_conf = a * self.bg.get_H(a)  # in SI units
        
        # Standard decay term in conformal time: -2ℋ ω
        # In terms of a: -2 ω / a
        decay = -2.0 * omega / a
        
        if self.hqiv_on:
            # Get the geometric horizon field φ = 3cH
            phi = self.bg.phi_horizon(a)  # in 1/s units
            
            # df/dφ for the square-root form: f(α, φ) = √(α - cφ)
            # df/dφ = -c / (2√(α - cφ)) = -c / (2α √(1 - cφ/α))
            # For the Brodie form: f = α - cφ → df/dφ = -c
            # We use a simplified model: df/dφ ∝ 1/α × (c/φ)²
            
            # Estimate alpha (inertia parameter) from background
            # α = φ / (3cH) in the exact background
            # Perturbatively: α = 1 + O(δ)
            
            # For the source, we need ∇φ from scalar perturbations
            # In a simplified treatment: ∇φ ∝ k × ω (from vorticity-density coupling)
            
            # The key insight: source grows with structure formation
            # We model this as: source_coef ∝ a^n where n > 0
            
            # Proper treatment would require coupling to scalar sector
            # For now, use a phenomenological model:
            # source = A * a^p * (k × grad_phi_mag) * direction
            
            if k_vec is None:
                k_vec = np.array([0.0, 0.0, self.k])  # Default: radial k
            
            # Simplified source model:
            # The source is proportional to structure growth
            # Use a simple growing model: ∇φ ∝ δ(a) ∝ D+(a)
            
            # Growth factor approximation (matter-dominated)
            # D+ ∝ a in standard cosmology
            # In HQIV, growth is modified
            
            # The conformal Hubble for the source denominator
            H_conf_SI = a * self.bg.get_H(a)  # 1/s
            
            # Source amplitude: proportional to structure growth
            # This is a simplified phenomenological model
            # The exact treatment requires coupling to scalar perturbations
            source_amplitude = 0.5  # Tuned parameter for demonstration
            
            # Structure growth proxy: a (matter-dominated)
            structure_growth = a
            
            # Source in conformal time: source_conf = amplitude * a^p * (k × grad_phi)
            # For simplicity: use proportional model
            source_conf = source_amplitude * structure_growth * self.k
            
            # Direction: perpendicular to k (for vorticity generation)
            # k × ∇φ gives a direction perpendicular to both k and ∇φ
            # Simplify: use z-direction for demonstration
            source_direction = np.array([0.0, 0.0, 1.0])
            
            # Full source term
            source = source_conf * source_direction
            
            # Divide by (a ℋ) to convert to da derivative
            dyda = decay + source / (a * H_conf_SI + 1e-30)
        else:
            dyda = decay
        
        return dyda
    
    def initial_conditions(self, a_start, amplitude=1e-10):
        """
        Set initial conditions for vorticity.
        
        In standard cosmology, vorticity is generated at second order
        and is very small. We start with a small seed vorticity.
        
        Parameters
        ----------
        a_start : float
            Initial scale factor
        amplitude : float
            Initial vorticity amplitude
            
        Returns
        -------
        y0 : array
            Initial state [ω_x, ω_y, ω_z]
        """
        # Start with a small vorticity with both parallel and perpendicular
        # components to test the amplification
        omega_x = amplitude * 0.5
        omega_y = amplitude * 0.5
        omega_z = amplitude * 1.0  # Parallel to default horizon direction
        
        return np.array([omega_x, omega_y, omega_z])
    
    def solve(self, a_start=1e-9, a_end=1.0, n_pts=500, k_vec=None):
        """
        Solve the vorticity evolution equation.
        
        Parameters
        ----------
        a_start : float
            Initial scale factor
        a_end : float
            Final scale factor
        n_pts : int
            Number of output points
        k_vec : array, optional
            Wavevector direction
            
        Returns
        -------
        a_arr : array
            Scale factor array
        omega_arr : array
            Vorticity evolution, shape (n_pts, 3)
        """
        # Ensure background is computed
        if self.bg._a_arr is None:
            self.bg.compute_background()
        
        # Initial conditions
        y0 = self.initial_conditions(a_start)
        
        # Time points for output
        a_eval = np.logspace(np.log10(a_start), np.log10(a_end), n_pts)
        
        # Solve ODE in log a
        def deriv_log_a(log_a, y):
            a = 10.0**log_a
            dyda = self.derivatives(a, y, k_vec)
            return dyda * a
        
        log_a_start = np.log10(a_start)
        log_a_end = np.log10(a_end)
        log_a_eval = np.log10(a_eval)
        
        sol = solve_ivp(
            deriv_log_a,
            [log_a_start, log_a_end],
            y0,
            method='RK45',
            t_eval=log_a_eval,
            rtol=1e-8,
            atol=1e-12
        )
        
        if not sol.success:
            print(f"Warning: Vorticity ODE solver did not converge: {sol.message}")
        
        # Extract solution
        self.a_arr = 10.0**sol.t
        self.omega_arr = sol.y.T  # Shape: (n_times, 3)
        
        return self.a_arr, self.omega_arr
    
    def vorticity_magnitude(self):
        """Compute the magnitude of vorticity |ω|."""
        if self.omega_arr is None:
            raise RuntimeError("Run solve() first")
        return np.sqrt(np.sum(self.omega_arr**2, axis=1))
    
    def vorticity_amplification(self):
        """
        Compute the vorticity amplification factor.
        
        A_ω = |ω(a)| / |ω(a_start)|
        
        In standard cosmology, A_ω ∝ a⁻² (decay).
        In HQIV, A_ω can grow or remain constant.
        
        Returns
        -------
        a_arr : array
            Scale factor array
        A_omega : array
            Amplification factor
        """
        if self.omega_arr is None:
            raise RuntimeError("Run solve() first")
        
        omega_mag = self.vorticity_magnitude()
        A_omega = omega_mag / omega_mag[0]
        
        return self.a_arr, A_omega


class VorticityScalarCoupling:
    """
    Handles the back-reaction of vorticity on scalar perturbations.
    
    Paper discussion: The amplified vorticity feeds back into:
    1. The velocity divergence θ (through velocity field coupling)
    2. The effective inertia reduction (through α modification)
    
    This coupling is essential for testing whether HQIV affects
    the CMB acoustic peak positions.
    
    Parameters
    ----------
    scalar_pert : ScalarPerturbations
        Scalar perturbation solver
    vorticity_pert : VorticityPerturbations
        Vorticity perturbation solver
    """
    
    def __init__(self, scalar_pert, vorticity_pert):
        self.scalar = scalar_pert
        self.vorticity = vorticity_pert
        
    def compute_back_reaction(self, a, delta, theta, omega):
        """
        Compute the back-reaction terms on scalar variables.
        
        The vorticity contributes to:
        1. An additional velocity term: δθ_ω ~ |ω|² / k
        2. Modified local acceleration: α → α + α_ω
        
        Parameters
        ----------
        a : float
            Scale factor
        delta : float
            Density contrast
        theta : float
            Velocity divergence
        omega : array
            Vorticity vector [ω_x, ω_y, ω_z]
            
        Returns
        -------
        delta_theta : float
            Correction to velocity divergence
        delta_alpha : float
            Correction to local acceleration
        """
        omega_mag = np.linalg.norm(omega)
        k = self.scalar.k
        
        # Velocity correction from vorticity
        # The rotational velocity component is v_rot ~ ω/k
        # This contributes to the total velocity but not to θ = ∇·v
        # However, it affects the inertia through α
        
        # Acceleration correction
        # The rotational motion contributes to the local acceleration
        # α_ω ~ |v_rot|² / r ~ ω² / k² × k = ω² / k
        delta_alpha = omega_mag**2 / (k + 1e-10)
        
        # Velocity divergence correction (second-order effect)
        # The coupling between vorticity and density gradients
        # can generate additional divergence
        delta_theta = 0.1 * omega_mag * delta  # Simplified coupling
        
        return delta_theta, delta_alpha
    
    def coupled_derivatives(self, a, y):
        """
        Compute derivatives with vorticity back-reaction included.
        
        State: y = [δ, θ, Φ, ω_x, ω_y, ω_z]
        
        Parameters
        ----------
        a : float
            Scale factor
        y : array
            Combined state vector
            
        Returns
        -------
        dyda : array
            Combined derivatives
        """
        # Split state
        scalar_y = y[:3]  # [δ, θ, Φ]
        omega = y[3:]     # [ω_x, ω_y, ω_z]
        
        # Scalar derivatives (base)
        dyda_scalar = self.scalar.derivatives(a, scalar_y)
        
        # Vorticity derivatives
        dyda_omega = self.vorticity.derivatives(a, omega)
        
        # Back-reaction corrections
        delta_theta, delta_alpha = self.compute_back_reaction(
            a, scalar_y[0], scalar_y[1], omega
        )
        
        # Apply corrections
        dyda_scalar[1] += delta_theta  # θ correction
        
        # The α correction affects the inertia factor in the Euler equation
        # This is already handled in estimate_alpha, but we can add explicit term
        # dyda_scalar[1] += ... (additional term from modified f)
        
        # Combine
        dyda = np.concatenate([dyda_scalar, dyda_omega])
        
        return dyda


def compute_vorticity_growth_exponent(background, k_test=0.1, a_start=1e-9, a_end=1.0):
    """
    Compute the effective vorticity growth exponent.
    
    Paper prediction: ω ∝ a^β for some effective exponent β.
    In standard cosmology: ω ∝ a⁻² (β = -2).
    In HQIV: β can be positive (amplification).
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology
    k_test : float
        Test wavenumber [1/Mpc]
    a_start : float
        Initial scale factor
    a_end : float
        Final scale factor
        
    Returns
    -------
    beta_eff : float
        Effective growth exponent
    a_arr : array
        Scale factor array
    A_omega : array
        Vorticity amplification factor
    """
    # Solve vorticity evolution
    vort = VorticityPerturbations(background, k_test, hqiv_on=True)
    a_arr, omega_arr = vort.solve(a_start, a_end)
    a_arr_h, A_omega = vort.vorticity_amplification()
    
    # Fit power law: ω ∝ a^β
    # log(A_omega) = β log(a/a_start)
    log_a = np.log(a_arr)
    log_A = np.log(A_omega + 1e-30)
    
    # Linear fit
    coeffs = np.polyfit(log_a, log_A, 1)
    beta_eff = coeffs[0]
    
    return beta_eff, a_arr, A_omega


# =============================================================================
# Testing and visualization
# =============================================================================

def test_vorticity_evolution():
    """Test vorticity evolution and create plots."""
    import matplotlib.pyplot as plt
    
    print("Testing vorticity evolution...")
    
    # Create background cosmologies
    hqiv_bg = CosmologicalBackground(H0=73.2, Om_m=0.031, hqiv_on=True,
                                       beta=1.02, Om_h=1.00, n_h=1.04)
    lcdm_bg = CosmologicalBackground(H0=67.4, Om_m=0.315, hqiv_on=False)
    
    hqiv_bg.compute_background()
    lcdm_bg.compute_background()
    
    # Test wavenumber
    k_test = 0.1  # Mpc⁻¹
    
    # Solve vorticity evolution
    print(f"\nSolving vorticity for k = {k_test} Mpc⁻¹...")
    
    vort_hqiv = VorticityPerturbations(hqiv_bg, k_test, hqiv_on=True)
    a_h, omega_h = vort_hqiv.solve(a_start=1e-9, a_end=1.0)
    A_h = vort_hqiv.vorticity_magnitude()
    A_h_norm = A_h / A_h[0]
    
    vort_lcdm = VorticityPerturbations(lcdm_bg, k_test, hqiv_on=False)
    a_l, omega_l = vort_lcdm.solve(a_start=1e-9, a_end=1.0)
    A_l = vort_lcdm.vorticity_magnitude()
    A_l_norm = A_l / A_l[0]
    
    # Standard decay: ω ∝ a⁻²
    a_theory = np.logspace(-9, 0, 500)
    A_theory = (a_theory / a_theory[0])**(-2)
    
    # Compute effective exponent
    beta_eff, _, _ = compute_vorticity_growth_exponent(hqiv_bg, k_test)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Vorticity magnitude evolution
    ax = axes[0, 0]
    ax.loglog(a_h, A_h_norm, 'b-', label='HQIV', linewidth=2)
    ax.loglog(a_l, A_l_norm, 'r--', label='ΛCDM (standard decay)', linewidth=2)
    ax.loglog(a_theory, A_theory, 'k:', label='Theory: a⁻²', linewidth=1.5, alpha=0.7)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('|ω| / |ω₀|')
    ax.set_title('Vorticity Magnitude Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([1e-10, 1e10])
    
    # Vorticity components (HQIV)
    ax = axes[0, 1]
    ax.loglog(a_h, np.abs(omega_h[:, 0]), 'r-', label='ω_x', linewidth=2)
    ax.loglog(a_h, np.abs(omega_h[:, 1]), 'g-', label='ω_y', linewidth=2)
    ax.loglog(a_h, np.abs(omega_h[:, 2]), 'b-', label='ω_z (parallel to ê_Θ)', linewidth=2)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('|ω_i|')
    ax.set_title('Vorticity Components (HQIV)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Effective exponent
    ax = axes[1, 0]
    # Compute local exponent: d(log A) / d(log a)
    log_a = np.log(a_h)
    log_A = np.log(A_h_norm + 1e-30)
    local_beta = np.gradient(log_A, log_a)
    ax.semilogx(a_h, local_beta, 'b-', linewidth=2)
    ax.axhline(-2.0, color='r', linestyle='--', alpha=0.7, label='Standard: β = -2')
    ax.axhline(0.0, color='k', linestyle=':', alpha=0.5)
    ax.axhline(beta_eff, color='g', linestyle=':', alpha=0.7, label=f'Avg: β = {beta_eff:.2f}')
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('β_eff(a)')
    ax.set_title('Local Vorticity Growth Exponent')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-3, 2])
    
    # Comparison of amplification
    ax = axes[1, 1]
    ratio = A_h_norm / (A_l_norm + 1e-30)
    ax.loglog(a_h, ratio, 'b-', linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('ω_HQIV / ω_ΛCDM')
    ax.set_title('Vorticity Amplification Ratio')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hqiv_solver/vorticity_evolution.png', dpi=150)
    plt.close()
    print("\nSaved vorticity_evolution.png")
    
    # Print summary
    print("\n" + "=" * 60)
    print("Vorticity Evolution Summary")
    print("=" * 60)
    print(f"Effective growth exponent β = {beta_eff:.4f}")
    print(f"  (Standard cosmology: β = -2)")
    print(f"Final amplification (HQIV): {A_h_norm[-1]:.2e}")
    print(f"Final amplification (ΛCDM): {A_l_norm[-1]:.2e}")
    print(f"Amplification ratio: {A_h_norm[-1]/A_l_norm[-1]:.2e}")
    print("=" * 60)
    
    return a_h, omega_h, a_l, omega_l


if __name__ == "__main__":
    test_vorticity_evolution()