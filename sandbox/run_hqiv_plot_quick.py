#!/usr/bin/env python3
"""Quick HQIV run (l_max=200) and plot multipole spectrum."""
import os
import sys

repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(repo, "CAMB"))

os.chdir(repo)
run_dir = os.path.join(repo, "camb_hqiv", "inifiles")
table_src = os.path.join(repo, "sandbox", "hqiv_Ha.txt")
if os.path.isfile(table_src):
    import shutil
    shutil.copy2(table_src, os.path.join(run_dir, "hqiv_Ha.txt"))
os.chdir(run_dir)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import camb

inifile = os.path.join(repo, "camb_hqiv", "inifiles", "params_hqiv.ini")
print("Loading HQIV params from:", inifile)
pars = camb.read_ini(inifile)
# Quick run
pars.set_for_lmax(200, lens_potential_accuracy=0)
print("Running CAMB (l_max=200)...")
results = camb.get_results(pars)

age = results.get_age()
print("Universe age today (Gyr): {:.4f}".format(age))

powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
cl_tt = powers["total"][:, 0]
ell = np.arange(2, len(cl_tt))
d_ell = ell * (ell + 1) * cl_tt[2:] / (2 * np.pi)

plt.figure(figsize=(10, 6))
plt.plot(ell, d_ell, color="C0", lw=1.2, label="HQIV CMB TT")
plt.xlim(2, ell[-1])
plt.ylim(0, max(8000, 1.05 * d_ell.max()))
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell(\ell+1)C_\ell/(2\pi)$  [$\mu$K$^2$]")
plt.title("HQIV CMB TT — tabulated H(a), no separate DE")
plt.legend()
plt.grid(True, alpha=0.3)
out = os.path.join(repo, "sandbox", "hqiv_cmb_multipole.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print("Saved:", out)
