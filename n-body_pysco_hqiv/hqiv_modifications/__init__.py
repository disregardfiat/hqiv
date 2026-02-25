"""
HQIV Modifications Package
==========================

This package implements the full action-derived HQIV physics for N-body simulations.
All equations are from paper/main.tex and paper/maxwell.tex.

Modules:
--------
- phi_field: Geometric horizon field φ(x) = 2c²/Θ_local(x)
- inertia_factor: Thermodynamic inertia reduction f(a_loc, φ)
- vorticity_source: Action-derived vorticity source term
- g_eff: Varying gravitational coupling G_eff(φ)
- bullet_ic: Bullet Cluster initial conditions
- background_wrapper: Interface to HQIV background solver

Key Equations (from paper):
---------------------------
1. Modified Einstein equation (paper/main.tex, Eq. after "Variation w.r.t. g^μν"):
   G_μν + γ(φ/c²) g_μν = (8π G_eff(φ)/c⁴) T_μν

2. Background Friedmann equation (paper/main.tex, Section 5):
   3H² - γH = 8π G_eff(H) (ρ_m + ρ_r)

3. Inertia factor (paper/main.tex, Eq. in Section 4.2):
   f(a_loc, φ) = max( a_loc / (a_loc + cφ/6), f_min )

4. Vorticity source (paper/main.tex, Eq. in Section 6):
   ∂ω/∂t + (v·∇)ω = (∂f/∂φ) (k × ∇φ) · ê_ω

5. Varying G (paper/main.tex, Eq. in Section 3):
   G(a) = G0 × (H(a)/H0)^α

Author: HQIV Team
"""

from .phi_field import (
    compute_phi_field,
    compute_theta_local,
    compute_Theta_local,
    phi_from_expansion_scalar,
)

from .inertia_factor import (
    inertia_reduction_factor,
    compute_local_acceleration,
    directional_inertia_factor,
)

from .vorticity_source import (
    vorticity_source_term,
    compute_vorticity_evolution,
    vorticity_growth_exponent,
)

from .g_eff import (
    effective_gravitational_constant,
    horizon_correction_term,
    modified_poisson_rhs,
)

from .bullet_ic import (
    BulletClusterIC,
    generate_nfw_halo,
    generate_collision_ic,
)

from .background_wrapper import (
    HQIVBackground,
    get_H_at_redshift,
    get_H_at_scale_factor,
    get_G_eff_at_redshift,
    get_default_background,
)

__all__ = [
    # phi_field
    'compute_phi_field',
    'compute_theta_local',
    'compute_Theta_local',
    'phi_from_expansion_scalar',
    # inertia_factor
    'inertia_reduction_factor',
    'compute_local_acceleration',
    'directional_inertia_factor',
    # vorticity_source
    'vorticity_source_term',
    'compute_vorticity_evolution',
    'vorticity_growth_exponent',
    # g_eff
    'effective_gravitational_constant',
    'horizon_correction_term',
    'modified_poisson_rhs',
    # bullet_ic
    'BulletClusterIC',
    'generate_nfw_halo',
    'generate_collision_ic',
    # background_wrapper
    'HQIVBackground',
    'get_H_at_redshift',
    'get_H_at_scale_factor',
    'get_G_eff_at_redshift',
    'get_default_background',
]

__version__ = '1.0.0'