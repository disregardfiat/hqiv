#!/usr/bin/env python3
"""
HQIV CAMB test script.
Run from CAMB repo root (or with PYTHONPATH including the built camb).
Uses params_hqiv.ini with HQIV = T and hqiv_Ha_file pointing to hqiv_Ha.txt.

Prints:
  - Universe age today (Gyr)
  - ℓ of first acoustic peak
  - Low-ℓ suppression (ℓ=2–30) in percent
"""
import os
import sys

def main():
    try:
        import camb
    except ImportError:
        print("CAMB not found. Build CAMB and install (e.g. pip install -e . from repo root).")
        sys.exit(1)

    # Paths: assume run from CAMB root or from camb_hqiv with repo as sibling
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    inifile = os.path.join(repo, "camb_hqiv", "inifiles", "params_hqiv.ini")
    if not os.path.isfile(inifile):
        inifile = os.path.join(repo, "inifiles", "params_hqiv.ini")
    if not os.path.isfile(inifile):
        inifile = "params_hqiv.ini"
    if not os.path.isfile(inifile):
        print("params_hqiv.ini not found. Copy camb_hqiv/inifiles/params_hqiv.ini to this dir or set path.")
        sys.exit(1)

    # HQIV table next to inifile or in current dir
    table = "hqiv_Ha.txt"
    if not os.path.isfile(table):
        table = os.path.join(repo, "sandbox", "hqiv_Ha.txt")
    if os.path.isfile(table):
        os.chdir(os.path.dirname(os.path.abspath(inifile)) or ".")
    else:
        print("Warning: hqiv_Ha.txt not found; ensure it is in the run directory or set hqiv_Ha_file in .ini")

    print("Loading params from:", inifile)
    pars = camb.read_ini(inifile)
    results = camb.get_results(pars)

    # Age today
    age = results.get_age()
    print("Universe age today (Gyr):", round(age, 4))

    # CMB Cl
    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
    cl_tt = powers["total"][:, 0]
    ell = list(range(len(cl_tt)))

    # First acoustic peak (rough: max in 180–280)
    start, end = 180, min(280, len(cl_tt) - 1)
    if end > start:
        peak_idx = start + max(range(end - start), key=lambda i: cl_tt[start + i])
        print("First acoustic peak: ℓ ≈", ell[peak_idx])
    else:
        print("First acoustic peak: (not enough l range)")

    # Low-ℓ suppression: compare to a reference (e.g. same spectrum without damping)
    low_ell = slice(2, 31)
    cl_low = cl_tt[low_ell]
    if cl_low.size and cl_low.mean() > 0:
        # If HQIV applies damping, mean is reduced; report relative drop
        ref = cl_tt[50:200].mean() if cl_tt[50:200].size else cl_low.mean()
        suppression = 100 * (1 - cl_low.mean() / ref) if ref > 0 else 0
        print("Low-ℓ (2–30) mean vs mid-ℓ reference: suppression ≈ {:.1f}%".format(suppression))
    else:
        print("Low-ℓ (2–30): (no data)")

    print("Done.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
