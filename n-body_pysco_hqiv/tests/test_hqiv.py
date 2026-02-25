#!/usr/bin/env python3
"""
HQIV N-body Module Tests
========================

Tests for the HQIV modifications package.

Run with:
    python3 tests/test_hqiv.py -v

Author: HQIV Team
"""

import sys
from pathlib import Path
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

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


def test_phi_field():
    """Test φ field computation."""
    print("\nTesting φ field...")
    
    # Create simple test case
    N = 32
    box_size = 50.0
    dx = box_size / N
    
    # Uniform density
    density = np.zeros((N, N, N))
    
    # Uniform velocity field (no expansion)
    velocity = np.zeros((3, N, N, N))
    
    # Compute φ
    H_a = 73.2 * 1000.0 / 3.0856775814913673e22  # H0 in 1/s
    a = 0.5
    
    phi, theta = compute_phi_field(density, velocity, H_a, a, dx)
    
    # Check that φ is computed
    assert phi is not None, "φ should not be None"
    assert phi.shape == (N, N, N), f"φ shape should be {(N, N, N)}, got {phi.shape}"
    
    # For uniform field, φ should be approximately constant
    phi_std = np.std(phi)
    assert phi_std < 1e-10, f"φ should be uniform for uniform input, std = {phi_std}"
    
    print(f"  φ range: [{phi.min():.2e}, {phi.max():.2e}]")
    print("  ✓ φ field test passed")
    
    return True


def test_inertia_factor():
    """Test inertia reduction factor."""
    print("\nTesting inertia factor...")
    
    # Test range of accelerations
    alpha = np.logspace(-14, -6, 100)
    
    # Typical φ value (c * H0 gives the horizon acceleration scale)
    H0 = 73.2 * 1000.0 / 3.0856775814913673e22
    phi = 2.99792458e8 * H0  # This is c * H0
    
    # Compute inertia factor
    f = inertia_reduction_factor(alpha, phi)
    
    # Check bounds
    assert np.all(f >= 0.01), "f should be >= f_min"
    assert np.all(f <= 1.0), "f should be <= 1"
    
    # Check monotonicity (higher acceleration → higher f)
    assert np.all(np.diff(f) >= 0), "f should increase with acceleration"
    
    # The minimum acceleration scale is a_min = chi * c * phi / 6
    # For phi = c * H0 ≈ 7e-10 m/s², a_min ≈ 0.172 * 3e8 * 7e-10 / 6 ≈ 0.006 m/s²
    # So we need accelerations >> 0.006 m/s² to get f → 1
    # And accelerations << 0.006 m/s² to get f → f_min
    
    f_high = inertia_reduction_factor(1.0, phi)  # 1 m/s² >> a_min
    f_low = inertia_reduction_factor(1e-14, phi)  # << a_min
    
    assert f_high > 0.9, f"f should be ~1 for very high acceleration, got {f_high}"
    assert f_low < 0.1, f"f should be small for low acceleration, got {f_low}"
    
    print(f"  f range: [{f.min():.4f}, {f.max():.4f}]")
    print(f"  f at high a (1.0 m/s²): {f_high:.4f}")
    print(f"  f at low a (1e-14 m/s²): {f_low:.4f}")
    print("  ✓ Inertia factor test passed")
    
    return True


def test_g_eff():
    """Test effective gravitational constant."""
    print("\nTesting G_eff...")
    
    H0 = 73.2 * 1000.0 / 3.0856775814913673e22
    
    # Test at different redshifts
    z = np.array([0, 1, 5, 10])
    a = 1.0 / (1 + z)
    
    # Approximate H(a) for matter-dominated
    H_a = H0 * np.sqrt(0.06 * a**(-3) + 0.94)
    
    # Compute G_eff
    G_eff = effective_gravitational_constant(H_a, H0, alpha=0.6)
    
    # G_eff should be >= 1 at high z (early times)
    assert G_eff[-1] > 1.0, f"G_eff at z=10 should be > 1, got {G_eff[-1]}"
    
    # G_eff should be ~1 at z=0
    assert 0.9 < G_eff[0] < 1.1, f"G_eff at z=0 should be ~1, got {G_eff[0]}"
    
    print(f"  G_eff at z=0: {G_eff[0]:.4f}")
    print(f"  G_eff at z=10: {G_eff[-1]:.4f}")
    print("  ✓ G_eff test passed")
    
    return True


def test_vorticity():
    """Test vorticity computation."""
    print("\nTesting vorticity...")
    
    N = 32
    box_size = 50.0
    dx = box_size / N
    
    # Create rotational velocity field
    x = np.linspace(0, box_size, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    center = box_size / 2
    
    velocity = np.zeros((3, N, N, N))
    omega_test = 0.1
    
    # Solid body rotation: v = ω × r
    # For rotation around z-axis: v_x = -ω(y-y0), v_y = ω(x-x0)
    velocity[0] = -omega_test * (Y - center)
    velocity[1] = omega_test * (X - center)
    
    # Compute vorticity
    from hqiv_modifications.vorticity_source import compute_vorticity_field
    omega = compute_vorticity_field(velocity, dx)
    
    # For solid body rotation, ω_z = 2 * omega_test
    # But FFT-based differentiation has numerical errors at boundaries
    # Check the mean in the central region
    omega_z_center = omega[2, N//4:3*N//4, N//4:3*N//4, N//4:3*N//4]
    omega_z_mean = np.mean(omega_z_center)
    expected = 2 * omega_test
    
    # Allow 50% tolerance due to FFT boundary effects
    assert np.abs(omega_z_mean - expected) < 0.5 * expected, \
        f"ω_z should be ~{expected}, got {omega_z_mean}"
    
    print(f"  Expected ω_z: {expected:.4f}")
    print(f"  Computed ω_z (center): {omega_z_mean:.4f}")
    print("  ✓ Vorticity test passed")
    
    return True


def test_bullet_ic():
    """Test Bullet Cluster initial conditions."""
    print("\nTesting Bullet Cluster IC...")
    
    import configparser
    
    # Small test configuration
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
            'npart_gas': '1000',
            'npart_galaxies': '500',
        }
    })
    
    # Generate IC
    ic = generate_collision_ic(test_config, seed=42)
    
    # Check output
    assert 'pos_gas' in ic, "IC should contain pos_gas"
    assert 'pos_gal' in ic, "IC should contain pos_gal"
    assert len(ic['pos_gas']) == 1000, f"Should have 1000 gas particles, got {len(ic['pos_gas'])}"
    assert len(ic['pos_gal']) == 500, f"Should have 500 galaxy particles, got {len(ic['pos_gal'])}"
    
    # Check positions are within box
    assert np.all(ic['pos_gas'] >= 0) and np.all(ic['pos_gas'] < 50.0), \
        "Gas positions should be within box"
    assert np.all(ic['pos_gal'] >= 0) and np.all(ic['pos_gal'] < 50.0), \
        "Galaxy positions should be within box"
    
    print(f"  Gas particles: {len(ic['pos_gas'])}")
    print(f"  Galaxy particles: {len(ic['pos_gal'])}")
    print(f"  Total mass: {ic['M_total']:.2e} Msun")
    print("  ✓ Bullet IC test passed")
    
    return True


def test_solver_integration():
    """Test integration of all solvers."""
    print("\nTesting solver integration...")
    
    # Create solvers
    phi_solver = PhiFieldSolver(H0=73.2, Omega_m=0.06, gamma=0.40,
                                 N_grid=32, box_size=50.0)
    inertia_solver = InertiaFactorSolver(H0=73.2, chi=0.172, f_min=0.01)
    gravity_solver = EffectiveGravitySolver(H0=73.2, Omega_m=0.06, gamma=0.40,
                                             alpha_G=0.6, N_grid=32, box_size=50.0)
    
    # Create test density
    density = np.zeros((32, 32, 32))
    density[16, 16, 16] = 10.0  # Point mass
    
    # Test gravity solver
    a = 0.5
    H_a = phi_solver.H_of_a(a)
    
    phi, acc, G_eff = gravity_solver.full_poisson_solve(density, a, H_a)
    
    assert phi is not None, "Potential should be computed"
    assert acc is not None, "Acceleration should be computed"
    assert G_eff > 0, "G_eff should be positive"
    
    print(f"  G_eff: {G_eff:.4f}")
    print(f"  Potential max: {phi.max():.4e}")
    print(f"  Acceleration max: {np.max(np.sqrt(np.sum(acc**2, axis=0))):.4e}")
    print("  ✓ Solver integration test passed")
    
    return True


def run_all_tests():
    """Run all tests."""
    print("="*60)
    print("HQIV N-body Module Tests")
    print("="*60)
    
    tests = [
        test_phi_field,
        test_inertia_factor,
        test_g_eff,
        test_vorticity,
        test_bullet_ic,
        test_solver_integration,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ Test failed: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"Results: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)