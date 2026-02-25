#!/usr/bin/env python3
"""
HQIV Background Wrapper — Quick Coherence Test
Run this to verify everything is working and producing sensible cosmology.
Uses gamma=0.4 (peak-alignment best-fit column) throughout.
"""

import sys
from pathlib import Path
import numpy as np

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("Note: matplotlib not available, skipping plots")

# Gamma from peak_alignment_scan3 best-fit (0.4 column)
GAMMA = 0.4
H0 = 73.2
OMEGA_M = 0.06

from hqiv_modifications.background_wrapper import (
    get_H_at_redshift,
    get_G_eff_at_redshift,
    HQIVBackground,
)

# ==================== TEST ====================
print("=== HQIV Background Coherence Test (γ = 0.4) ===\n")

# 1. Universe age (should be ~32–34 Gyr) — integrated from H(z)
print("Integrating universe age from z=∞ to z=0...")
# Proper time integral: t = ∫ da / (a * H(a))
# from a = 1e-4 (early universe) to a = 1 (today)
a = np.logspace(-4, 0, 1000)  # a goes from small to large
H = get_H_at_redshift(1 / a - 1, gamma=GAMMA, H0=H0, Omega_m=OMEGA_M)
H_SI = H * 1000.0 / 3.0856775814913673e22  # Convert to 1/s

# dt = da / (a * H)
dt = np.diff(a) / (a[:-1] * H_SI[:-1])
age_integrated_Gyr = np.sum(dt) / (3600 * 24 * 365.25 * 1e9)
print(f"  Age (integrated): {age_integrated_Gyr:.4f} Gyr   ← ~19 Gyr for γ=0.4, H0=73.2, Ω_m=0.06")

# 2. Age from wrapper / hqiv_solver (must match)
bg = HQIVBackground(H0=H0, Omega_m=OMEGA_M, gamma=GAMMA)
age_wrapper_Gyr = bg.age_today()
print(f"  Age (wrapper):    {age_wrapper_Gyr:.4f} Gyr")

# Verify ages match (within 1% or 0.3 Gyr)
age_diff = abs(age_integrated_Gyr - age_wrapper_Gyr)
age_ok = age_diff < max(0.3, 0.01 * age_wrapper_Gyr)
if age_ok:
    print(f"  → Ages match (Δ = {age_diff:.4f} Gyr).\n")
else:
    print(f"  → WARNING: age mismatch Δ = {age_diff:.4f} Gyr\n")
    sys.exit(1)

# 3. Key epochs
print("Key values (γ=0.4, H0=73.2, Ω_m=0.06):")
for z in [0, 0.296, 1, 10, 1090]:
    H = get_H_at_redshift(z, gamma=GAMMA, H0=H0, Omega_m=OMEGA_M)
    G = get_G_eff_at_redshift(z, gamma=GAMMA, H0=H0, Omega_m=OMEGA_M)
    print(f"  z = {z:6.1f} | H = {H:7.2f} km/s/Mpc | G_eff = {G:6.3f} G₀")

# 4. Plots (if matplotlib available)
if HAS_MPL:
    z_plot = np.logspace(-2, 4, 500)
    H_plot = get_H_at_redshift(z_plot, gamma=GAMMA, H0=H0, Omega_m=OMEGA_M)
    G_plot = get_G_eff_at_redshift(z_plot, gamma=GAMMA, H0=H0, Omega_m=OMEGA_M)

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax[0].loglog(z_plot, H_plot, 'b-', lw=2, label='HQIV')
    ax[0].set_ylabel(r'H(z) [km s⁻¹ Mpc⁻¹]')
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    ax[1].semilogx(z_plot, G_plot, 'r-', lw=2, label='G$_{\\rm eff}$(z)')
    ax[1].set_xlabel('Redshift z')
    ax[1].set_ylabel(r'G$_{\rm eff}$ / G₀')
    ax[1].grid(True, alpha=0.3)
    ax[1].legend()

    plt.suptitle('HQIV Background Evolution (γ=0.4, H0=73.2, Ω_m=0.06)', fontsize=14)
    plt.tight_layout()
    plt.savefig('hqiv_background_coherence_test.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\n✅ Test complete!")
    print("   → Plot saved as hqiv_background_coherence_test.png")
else:
    print("\n✅ Test complete!")
    print("   → Install matplotlib to see plots")
    print("   → Age and key values printed above")
    print("   → Coherent if integrated and wrapper ages match (~19 Gyr for this cosmology).")