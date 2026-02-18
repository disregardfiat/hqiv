#!/usr/bin/env python3
"""
Short Python wrapper: load (modified) CAMB, run with HQIV params, produce TT spectrum and save plot.

Usage:
  python hqiv_camb_wrapper.py [path_to_params.ini]

If no .ini given, uses camb_hqiv/inifiles/params_hqiv.ini (with HQIV = T and hqiv_Ha_file).
Expects hqiv_Ha.txt in the current directory or path set in the .ini.
"""
import os
import sys

def main():
    try:
        import camb
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        print("Need: camb, numpy, matplotlib. Error:", e)
        sys.exit(1)

    base = os.path.dirname(os.path.abspath(__file__))
    inifile = sys.argv[1] if len(sys.argv) > 1 else None
    if not inifile:
        inifile = os.path.join(base, "inifiles", "params_hqiv.ini")
    if not os.path.isfile(inifile):
        inifile = "params_hqiv.ini"
    if not os.path.isfile(inifile):
        print("No params .ini found. Pass path or copy params_hqiv.ini here.")
        sys.exit(1)

    pars = camb.read_ini(inifile)
    results = camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
    cl_tt = powers["total"][:, 0]
    ell = np.arange(len(cl_tt))

    outfig = os.path.join(os.getcwd(), "hqiv_cmb_tt.png")
    plt.figure(figsize=(10, 6))
    plt.plot(ell, ell * (ell + 1) * cl_tt / (2 * np.pi), label="HQIV TT", lw=1.5)
    plt.xlim(2, 1200)
    plt.ylim(0, None)
    plt.xlabel(r"Multipole $\ell$")
    plt.ylabel(r"$\ell(\ell+1)C_\ell/(2\pi)$ [$\mu$K$^2$]")
    plt.title("HQIV CMB TT spectrum")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(outfig, dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved:", outfig)
    return 0

if __name__ == "__main__":
    sys.exit(main())
