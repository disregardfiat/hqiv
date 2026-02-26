#!/usr/bin/env python3
"""
Plot (1) H(χ) at "now" — first-order radial expansion (paper: H(χ) = H_loc − |∂H/∂χ| χ + O(χ²);
     gradient ∝ sectional curvature K ≈ −H_loc² of H³ fibre, damped by γ ≈ 0.40);
     (2) Apparent (experienced) lookback time at z and volume-averaged H(z)/H₀.

Paper: observer-centric profile; background and CLASS use only ⟨φ⟩ = c⟨H(a)⟩. Fiducial:
H_loc = 73 km/s/Mpc, H_vol = 16.1 km/s/Mpc (⟨H⟩ at z=0), conformal age 49597 Mpc.
Output: paper/radial_H_and_age.pdf
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Fiducial (from CLASS paper_run full_run.log)
H_INPUT_MPC = 2.435018e-04    # 1/Mpc -> 73 km/s/Mpc
H_ACTUAL_MPC = 5.365041e-05   # 1/Mpc -> 16.1 km/s/Mpc
CONFORMAL_AGE_MPC = 49597.354274
AGE_GYR = 51.175906           # wall-clock age
COMPRESSION = 3.96            # wall-clock / apparent ΛCDM age at z=0
C_KM_S = 2.998e5              # c in km/s

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PAPER_RUN = os.path.join(SCRIPT_DIR, 'paper_run')
BG_FILE = os.path.join(PAPER_RUN, 'hqiv_min00_background.dat')
OUT_DIR = os.path.join(SCRIPT_DIR, '..', 'paper')
OUT_PDF = os.path.join(OUT_DIR, 'radial_H_and_age.pdf')


def main():
    # --- Panel 1: H(χ) at "now" (paper's first-order gradient) ---
    # Paper: H(χ) = H_loc − |∂H/∂χ|_{χ=0} χ + O(χ²); |∂H/∂χ| ∝ K ≈ −H_loc² (damped by γ).
    # Fiducial normalization: set gradient so H(χ_horizon) = H_vol (⟨H⟩); illustrative.
    chi_horizon = CONFORMAL_AGE_MPC
    H_loc = H_INPUT_MPC * C_KM_S   # km/s/Mpc (observer, χ=0)
    H_vol = H_ACTUAL_MPC * C_KM_S  # km/s/Mpc (volume-averaged; CLASS H(z=0))
    chi_mpc = np.linspace(0, chi_horizon, 300)
    # First-order profile: H(χ_horizon) ≈ H_vol
    dH_dchi = (H_loc - H_vol) / chi_horizon
    H_chi = H_loc - dH_dchi * chi_mpc

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

    # Left: H(chi) and H(chi)/H_loc
    ax1.plot(chi_mpc, H_chi, 'b-', lw=1.5, label=r'$H(\chi)$ [km/s/Mpc]')
    ax1.axhline(H_vol, color='gray', ls='--', alpha=0.7, label=r'$H_{\rm vol}$ (16.1)')
    ax1.axhline(H_loc, color='gray', ls=':', alpha=0.7, label=r'$H_{\rm loc}$ (73)')
    ax1.set_xlabel(r'Comoving radial distance $\chi$ [Mpc]')
    ax1.set_ylabel(r'$H(\chi)$ [km/s/Mpc]')
    ax1.set_title(r'Radial $H$ at "now" (observer $\chi=0$)')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, None)

    # Same panel, ratio H(chi)/H_loc on secondary y-axis
    ax1_twin = ax1.twinx()
    ax1_twin.plot(chi_mpc, H_chi / H_loc, 'darkgreen', lw=1, alpha=0.8, label=r'$H(\chi)/H_{\rm loc}$')
    ax1_twin.set_ylabel(r'$H(\chi)/H_{\rm loc}$', color='darkgreen')
    ax1_twin.tick_params(axis='y', labelcolor='darkgreen')
    ax1_twin.set_ylim(0, 1.05)

    # --- Panel 2: Apparent (experienced) age at z, and H(z)/H0 ---
    if not os.path.isfile(BG_FILE):
        ax2.text(0.5, 0.5, 'No background file', ha='center', va='center', transform=ax2.transAxes)
    else:
        bg = np.loadtxt(BG_FILE)
        if bg.ndim == 1:
            bg = bg.reshape(1, -1)
        z = bg[:, 0]
        proper_time_gyr = bg[:, 1]   # cosmic time at that z (wall-clock)
        H_z_mpc = bg[:, 3]           # H(z) in 1/Mpc
        age_today = proper_time_gyr[z <= 1e-6][-1] if np.any(z <= 1e-6) else AGE_GYR
        # Wall-clock lookback (HQIV global time)
        lookback_gyr = age_today - proper_time_gyr
        # Apparent (experienced) lookback: compressed by factor ~3.96
        lookback_app_gyr = lookback_gyr / COMPRESSION
        H0_mpc = H_ACTUAL_MPC
        H_over_H0 = H_z_mpc / H0_mpc
        # Log-log panel from z ~ 1e-2 up to recombination (~1100)
        mask = (z > 1e-2) & (z <= 1.2e3) & (lookback_app_gyr > 0)
        z_plot = z[mask]
        lookback_app_plot = lookback_app_gyr[mask]
        H_ratio_plot = H_over_H0[mask]
        # Apparent (experienced) lookback
        ax2.plot(z_plot, lookback_app_plot, 'b-', lw=1.5, label='Lookback (apparent)')
        # Horizontal guide: ~13.8 Gyr apparent age
        ax2.axhline(13.8, color='gray', ls='--', alpha=0.7, label=r'$\sim$13.8 Gyr (ΛCDM chronometers)')
        ax2.set_xlabel(r'Redshift $z$')
        ax2.set_ylabel(r'Apparent lookback time [Gyr]', color='b')
        ax2.tick_params(axis='y', labelcolor='b')
        ax2.set_title(r'Experienced age at $z$ and $H(z)/H_0$')
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlim(1e-2, 1.2e3)
        ax2.set_ylim(1e-2, 2e2)
        ax2_twin = ax2.twinx()
        ax2_twin.plot(z_plot, H_ratio_plot, 'darkgreen', lw=1, alpha=0.8, label=r'$H(z)/H_0$')
        ax2_twin.set_ylabel(r'$H(z)/H_0$ (volume-avg.)', color='darkgreen')
        ax2_twin.tick_params(axis='y', labelcolor='darkgreen')
        ax2_twin.set_ylim(0, None)

    fig.suptitle(r'HQIV fiducial: radial $H(\chi)$ and lookback time', fontsize=10)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close()
    print('Saved', OUT_PDF)


if __name__ == '__main__':
    main()
