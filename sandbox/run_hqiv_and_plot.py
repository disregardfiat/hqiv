#!/usr/bin/env python3
"""Run HQIV CAMB and plot CMB TT multipole spectrum (D_ell vs ell)."""
import os
import sys

repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(repo, "CAMB"))

os.chdir(repo)
inifile = os.path.join(repo, "camb_hqiv", "inifiles", "params_hqiv.ini")
table = os.path.join(repo, "sandbox", "hqiv_Ha.txt")
run_dir = os.path.join(repo, "camb_hqiv", "inifiles")
if os.path.isfile(table):
    # copy table so CAMB finds it when run from inifiles dir
    import shutil
    dest = os.path.join(run_dir, "hqiv_Ha.txt")
    if os.path.dirname(os.path.abspath(inifile)) != os.path.dirname(os.path.abspath(table)):
        shutil.copy2(table, dest)
os.chdir(run_dir)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import camb

print("Loading HQIV params from:", inifile)
pars = camb.read_ini(inifile)
print("Running CAMB (HQIV background)...")
results = camb.get_results(pars)

age = results.get_age()
print("Universe age today (Gyr, HQIV prediction): {:.4f}".format(age))

powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
cl_tt = powers["total"][:, 0]
ell = np.arange(2, len(cl_tt))
cl_tt = cl_tt[2:]
# D_ell = ell(ell+1) C_ell / (2 pi)
d_ell = ell * (ell + 1) * cl_tt / (2 * np.pi)

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
