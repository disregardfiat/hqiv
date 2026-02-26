#!/usr/bin/env python3
"""
HQIV UNITY EXPLORER — Everything in pure [0,1] geometry
Energy 0 → 1 (Planck), horizon = 1/E', angle of time rotates Fano triples
"""

import numpy as np

# =============================================================================
# PURE UNITY CONSTANTS
# =============================================================================
E_PL = 1.0                    # everything normalized to Planck = 1
T_QCD_PRIME = 1.8 / 1.22e19   # ~1.475e-19  (QCD in Planck units)

FANO_LINES = [(1,2,4), (2,5,7), (3,6,1), (4,7,2), (5,1,6), (6,3,4), (7,4,5)]

def unity_energy(E_GeV):
    """Normalize energy to [0,1]"""
    return E_GeV / 1.22e19

def horizon_prime(E_prime):
    """Θ' = 1 / E' — the horizon built by the moment of the energy"""
    return 1.0 / np.maximum(E_prime, 1e-300)

def phi_prime(E_prime):
    """ϕ' = 2 / Θ' = 2 E'"""
    return 2.0 * E_prime

def angle_of_time(E_prime):
    """The angle of time, pure geometry"""
    return np.arctan(E_prime) * (np.pi / 2)

def curvature_imprint_prime(m_prime):
    """Normalized curvature imprint (same mechanism as Ω_k^true)"""
    return 0.0098 / (m_prime + 1) * (1 + 0.6 * np.log(1.0 / (m_prime + 1))) * 4.84e5

# =============================================================================
# RUN THE EXPLORER
# =============================================================================
print("HQIV UNITY EXPLORER — All in pure geometry [0,1]")
print("=" * 75)

for E_GeV in [1.8, 5, 30, 100, 1000, 1.22e16]:
    E_p = unity_energy(E_GeV)
    Θ_p = horizon_prime(E_p)
    ϕ_p = phi_prime(E_p)
    θ = angle_of_time(E_p)
    
    print(f"\nEnergy = {E_GeV:8.1e} GeV → E' = {E_p:.2e}")
    print(f"  Horizon Θ'     = {Θ_p:.2e}")
    print(f"  ϕ'             = {ϕ_p:.2e}")
    print(f"  Angle of time  = {θ:.4f} rad  ({θ*180/np.pi:.1f}°)")
    
    shift = int(θ / (np.pi / 7)) % 7
    print(f"  Fano rotation  = +{shift} steps")
    print("  Active triples : ", end='')
    for line in FANO_LINES:
        rotated = tuple((x + shift - 1) % 7 + 1 for x in line)
        print(rotated, end=' ')
    print()