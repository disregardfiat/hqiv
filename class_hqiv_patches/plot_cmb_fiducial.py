#!/usr/bin/env python3
"""
Plot CMB TT power spectrum D_ell = ell(ell+1)C_ell/(2pi) vs multipole ell
at the fiducial cost-minimum parameters. Saves figure for paper.
CLASS output is dimensionless [l(l+1)/2pi] C_l; we plot vs l (2 to l_max).
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Conformal-time compression (wall-clock / apparent age); used only for top axis
TIME_COMPRESSION = 3.96
# T_cmb in muK for Planck D_ell (muK^2) -> dimensionless
T_CMB_MUK = 2.725e6
# Set True to show top axis with observer-centric multipole ell/compression (illustrative only)
SHOW_TOP_AXIS = True

# Paths: paper_run has the fiducial full run; fallback to run_at_minimum
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PAPER_RUN = os.path.join(SCRIPT_DIR, 'paper_run')
RUN_AT_MIN = os.path.join(SCRIPT_DIR, 'run_at_minimum')
RUN_DIR = PAPER_RUN if os.path.isdir(PAPER_RUN) and os.path.isfile(os.path.join(PAPER_RUN, 'hqiv_min00_cl.dat')) else RUN_AT_MIN
CL_FILE = os.path.join(RUN_DIR, 'hqiv_min00_cl.dat')
PLANCK_BINNED_FILE = os.path.join(SCRIPT_DIR, 'planck2018_tt_binned.dat')
PLANCK_LOWELL_FILE = os.path.join(SCRIPT_DIR, 'planck2018_tt_lowell.dat')
OUT_DIR = os.path.join(SCRIPT_DIR, '..', 'paper')
OUT_PDF = os.path.join(OUT_DIR, 'cmb_spectrum_fiducial.pdf')

def main():
    if not os.path.isfile(CL_FILE):
        raise SystemExit('CL file not found: ' + CL_FILE)
    data = np.loadtxt(CL_FILE, comments='#')
    if data.ndim == 1:
        data = data.reshape(-1, 2)
    ell = data[:, 0].astype(int)
    # Column 2 is dimensionless l(l+1)C_l/2pi (CLASS header)
    d_ell = data[:, 1]

    # Planck 2018 TT: low-ell (axis of evil) and binned high-ell; D_ell muK^2 -> dimensionless
    def load_planck_dl(path):
        if not os.path.isfile(path):
            return None, None, None
        pl = np.loadtxt(path, comments='#')
        if pl.ndim == 1:
            pl = pl.reshape(-1, 4)
        if pl.ndim != 2 or pl.shape[1] < 2:
            return None, None, None
        ell = pl[:, 0]
        d_muk2 = pl[:, 1]
        d_dim = d_muk2 / (T_CMB_MUK ** 2)
        err = (pl[:, 2] + pl[:, 3]) / (2.0 * T_CMB_MUK ** 2) if pl.shape[1] >= 4 else None
        return ell, d_dim, err

    planck_lo_ell, planck_lo_d, planck_lo_err = load_planck_dl(PLANCK_LOWELL_FILE)
    planck_hi_ell, planck_hi_d, planck_hi_err = load_planck_dl(PLANCK_BINNED_FILE)

    fig, ax = plt.subplots(1, 1, figsize=(7.0, 4.0))
    # Axis-of-evil region (low ell): shaded band + Planck low-ell points
    ax.axvspan(2, 30, alpha=0.12, color='green', zorder=0)
    if planck_lo_ell is not None:
        if planck_lo_err is not None:
            ax.errorbar(planck_lo_ell, planck_lo_d, yerr=planck_lo_err, fmt='o', color='#2ca02c',
                       markersize=4, capsize=1.5, alpha=0.9, label='Planck 2018 TT (low-$\\ell$, Commander)', zorder=3)
        else:
            ax.plot(planck_lo_ell, planck_lo_d, 'o', color='#2ca02c', markersize=4, alpha=0.9,
                   label='Planck 2018 TT (low-$\\ell$)', zorder=3)
    # Planck high-ell (binned)
    if planck_hi_ell is not None:
        if planck_hi_err is not None:
            ax.errorbar(planck_hi_ell, planck_hi_d, yerr=planck_hi_err, fmt='.', color='#2ca02c',
                       markersize=3, capsize=0, alpha=0.7, label='Planck 2018 TT (binned)', zorder=2)
        else:
            ax.plot(planck_hi_ell, planck_hi_d, '.', color='#2ca02c', markersize=3, alpha=0.8,
                   label='Planck 2018 TT (binned)', zorder=2)
    # Raw HQIV spectrum (volume-averaged background; peaks at higher ell than Planck)
    ax.plot(ell, d_ell, color='#1f77b4', linewidth=1.2, label=r'HQIV fiducial (volume-avg.)')
    ax.set_xscale('linear')
    ax.set_xlabel(r'Multipole $\ell$ (volume-averaged)')
    ax.set_ylabel(r'$D_\ell \equiv \ell(\ell+1)C_\ell/(2\pi)$ [dimensionless]')
    ax.set_xlim(2, ell.max())
    ax.set_ylim(0, None)
    ymax = ax.get_ylim()[1]
    ax.text(8, ymax * 0.92, r'axis of evil ($\ell \lesssim 20$)', fontsize=8, color='darkgreen', alpha=0.9)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, alpha=0.3)
    ax.set_title(r'CMB TT power spectrum: Planck 2018 vs HQIV fiducial')
    # Optional top axis: observer-centric multipole (illustrative; full correction in perturbation code)
    if SHOW_TOP_AXIS:
        ax_top = ax.secondary_xaxis('top', functions=(lambda l: l / TIME_COMPRESSION, lambda l_eff: l_eff * TIME_COMPRESSION))
        ax_top.set_xlabel(r'$\ell_{\rm eff}$ (observer-centric, $\ell/3.96$)', fontsize=9)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close()
    print('Saved', OUT_PDF)

if __name__ == '__main__':
    main()
