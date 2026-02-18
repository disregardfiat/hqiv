#!/usr/bin/env python3
"""
Full covariant HQIV (McCulloch Quantised Inertia).

Uses params with HQIV = T, HQIV_covariant = T, hqiv_beta = 1.02, and hqiv_Ha.txt
generated with dynamic horizon term: A_eff = A_std + beta*H(t)^2 (target age ~17–20 Gyr).

Prints:
  - Universe age today (Gyr) — must be ~17–20 Gyr
  - First acoustic peak ℓ
  - Low-ℓ suppression (ℓ=2–30) %
  - σ8
  - D(z=14) / D_ΛCDM

Requires: built CAMB with HQIV + HQIV_covariant; hqiv_Ha.txt in run dir or path in .ini.
"""
from __future__ import division

import os
import sys

def main():
    try:
        import camb
        import numpy as np
    except ImportError:
        print("CAMB/numpy not found. Build CAMB with HQIV covariant patches and install.")
        sys.exit(1)

    base = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(base)
    inifile = os.path.join(base, "inifiles", "params_hqiv_covariant.ini")
    if not os.path.isfile(inifile):
        inifile = os.path.join(base, "inifiles", "params_hqiv.ini")
    if not os.path.isfile(inifile):
        inifile = os.path.join(repo, "inifiles", "params_hqiv.ini")
    if not os.path.isfile(inifile):
        inifile = "params_hqiv.ini"
    if not os.path.isfile(inifile):
        print("No params .ini found. Copy params_hqiv_covariant.ini or params_hqiv.ini.")
        sys.exit(1)

    table = "hqiv_Ha.txt"
    if not os.path.isfile(table):
        table = os.path.join(repo, "sandbox", "hqiv_Ha.txt")
    if os.path.isfile(table):
        run_dir = os.path.dirname(os.path.abspath(inifile))
        if run_dir:
            os.chdir(run_dir)
    else:
        print("Warning: hqiv_Ha.txt not found; set hqiv_Ha_file in .ini.")

    print("Loading HQIV (covariant) params from:", inifile)
    print("HQIV: Omega_total=1 from horizon equilibrium (no separate DE); age is model prediction.")
    pars = camb.read_ini(inifile)
    pars.WantTransfer = True
    pars.Transfer.PK_num_redshifts = 2
    pars.Transfer.PK_redshifts = [14.0, 0.0]

    results = camb.get_results(pars)

    age = results.get_age()
    print("Universe age today (Gyr, HQIV prediction): {:.4f}".format(age))

    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
    cl_tt = powers["total"][:, 0]
    ell = np.arange(len(cl_tt))

    start, end = 180, min(280, len(cl_tt) - 1)
    if end > start:
        peak_idx = start + np.argmax(cl_tt[start:end])
        print("First acoustic peak: ℓ ≈ {:d}".format(ell[peak_idx]))
    else:
        print("First acoustic peak: (not enough l range)")

    low_ell = slice(2, 31)
    cl_low = cl_tt[low_ell]
    if cl_low.size and cl_low.mean() > 0:
        ref = cl_tt[50:200].mean() if cl_tt[50:200].size else cl_low.mean()
        suppression = 100 * (1 - cl_low.mean() / ref) if ref > 0 else 0
        print("Low-ℓ (2–30) suppression: {:.1f}%".format(suppression))
    else:
        print("Low-ℓ (2–30): (no data)")

    try:
        sigma8 = results.get_sigma8_0()
        print("σ8 today: {:.4f}".format(sigma8))
    except Exception as e:
        print("σ8 today: (failed – ensure WantTransfer and z=0 in PK_redshifts)", e)

    k_small = 1e-3
    z_ev = [0.0, 14.0]
    try:
        ev = results.get_redshift_evolution(k_small, z_ev, vars=["delta_tot"])
        if ev.shape[0] == 1:
            ev = ev[0]
        iz0 = 0
        iz14 = 1
        if z_ev[0] > z_ev[1]:
            iz0, iz14 = 1, 0
        d0 = float(ev[iz0])
        d14 = float(ev[iz14])
        if abs(d0) > 1e-30:
            D14_HQIV = d14 / d0
        else:
            D14_HQIV = np.nan
    except Exception as e:
        D14_HQIV = np.nan
        print("D(z=14) HQIV: (failed)", e)

    print("D(z=14) HQIV (delta(14)/delta(0)): {:.6f}".format(D14_HQIV))

    print("Running ΛCDM for D(z=14) comparison...")
    pars_lcdm = camb.read_ini(inifile)
    pars_lcdm.HQIV = False
    if hasattr(pars_lcdm, "HQIV_covariant"):
        pars_lcdm.HQIV_covariant = False
    pars_lcdm.omch2 = 0.12
    pars_lcdm.ombh2 = 0.048
    pars_lcdm.H0 = pars.H0
    pars_lcdm.WantTransfer = True
    pars_lcdm.Transfer.PK_num_redshifts = 2
    pars_lcdm.Transfer.PK_redshifts = [14.0, 0.0]

    results_lcdm = camb.get_results(pars_lcdm)
    try:
        ev_lcdm = results_lcdm.get_redshift_evolution(k_small, z_ev, vars=["delta_tot"])
        if ev_lcdm.shape[0] == 1:
            ev_lcdm = ev_lcdm[0]
        d0_lcdm = float(ev_lcdm[iz0])
        d14_lcdm = float(ev_lcdm[iz14])
        if abs(d0_lcdm) > 1e-30:
            D14_LCDM = d14_lcdm / d0_lcdm
        else:
            D14_LCDM = np.nan
    except Exception as e:
        D14_LCDM = np.nan
        print("D(z=14) ΛCDM: (failed)", e)

    print("D(z=14) ΛCDM: {:.6f}".format(D14_LCDM))
    if not (np.isnan(D14_HQIV) or np.isnan(D14_LCDM)) and D14_LCDM != 0:
        ratio = D14_HQIV / D14_LCDM
        print("D(z=14) HQIV / D(z=14) ΛCDM: {:.6f}".format(ratio))
    print("Done.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
