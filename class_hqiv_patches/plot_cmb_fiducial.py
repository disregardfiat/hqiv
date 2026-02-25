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

# Paths: run_at_minimum has the fiducial cl file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.path.join(SCRIPT_DIR, 'run_at_minimum')
CL_FILE = os.path.join(RUN_DIR, 'hqiv_min00_cl.dat')
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

    fig, ax = plt.subplots(1, 1, figsize=(6, 3.5))
    ax.plot(ell, d_ell, color='#1f77b4', linewidth=1.2, label=r'HQIV fiducial ($\Omega_m=0.0191$, cost min.)')
    ax.set_xscale('log')
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$D_\ell \equiv \ell(\ell+1)C_\ell/(2\pi)$ [dimensionless]')
    ax.set_xlim(2, ell.max())
    ax.set_ylim(0, None)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, alpha=0.3)
    ax.set_title(r'CMB TT power spectrum at fiducial parameters')
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close()
    print('Saved', OUT_PDF)

if __name__ == '__main__':
    main()
