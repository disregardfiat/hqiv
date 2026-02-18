"""Quick HQIV test with clean table: small l_max first to confirm no DVERK."""
import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isfile("hqiv_Ha.txt"):
    print("Run generate_hqiv_background.py first.")
    sys.exit(1)

import camb

# Load HQIV params from ini (must have HQIV=T and hqiv_Ha_file)
inifile = os.path.join(os.path.dirname(__file__), "..", "camb_hqiv", "inifiles", "params_hqiv.ini")
if os.path.isfile(inifile):
    pars = camb.read_ini(inifile)
else:
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=70, ombh2=0.048 * (70/100)**2, omch2=0, omk=0)
    pars.set_for_lmax(200, lens_potential_accuracy=0)
    # If your camb has HQIV on the params struct, set them here when not using ini
    if hasattr(pars, "HQIV"):
        pars.HQIV = True
        pars.hqiv_Ha_file = "hqiv_Ha.txt".ljust(1024)[:1024]

# Small test first
pars.set_for_lmax(200, lens_potential_accuracy=0)

print("Running CAMB (l_max=200) with clean hqiv_Ha.txt ...")
results = camb.get_results(pars)
print("  get_results succeeded.")
print("Universe age today (Gyr):", round(results.get_age(), 4))
powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
cl_tt = powers["total"][:, 0]
ell = list(range(len(cl_tt)))
start, end = 180, min(280, len(cl_tt) - 1)
if end > start:
    peak_idx = start + max(range(end - start), key=lambda i: cl_tt[start + i])
    print("First acoustic peak: ℓ ≈", ell[peak_idx])
print("✅ Safe test passed — ready for full l_max=2500 + horizon cutoff.")
