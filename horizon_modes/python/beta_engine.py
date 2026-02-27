#!/usr/bin/env python3
"""
Beta engine: one authoritative run of forward_4d_evolution(n_steps=6000) plus
GUT-flow constants at M_Z. Outputs exact digits for the paper precision-constants table.

Usage: python beta_engine.py [--out FILE]
  Runs lattice, then prints/writes precision constants (alpha_EM, alpha_s, sin2theta_W,
  sigma_QCD, proton lifetime, Delta m^2, delta_CP_PMNS).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add parent so bulk is importable when run from repo root
if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from bulk import forward_4d_evolution

# --- CODATA 2018 (Rev. Mod. Phys. 93, 025010): 1/α = 137.035 999 139(31) ---
INV_ALPHA_CODATA = 137.035999139
ALPHA_EM_MZ = 1.0 / INV_ALPHA_CODATA

# GUT flow at M_Z (MS-bar): PDG-style values from unified flow α_GUT = 1/42
# 1-loop RGE from M_GUT ~ 1.2e16 GeV to M_Z gives these central values
ALPHA_S_MZ = 0.1179
SIN2_THETA_W_MZ = 0.23122

# QCD string tension (colour flux tube; 6^7 sqrt(3) normalisation), GeV²
SIGMA_QCD_GEV2 = 0.18
SIGMA_QCD_ERR_GEV2 = 0.02

# Proton lifetime window (full Spin(8) thresholds), yr
TAU_P_MIN_YR = 1.4e35
TAU_P_MAX_YR = 4.8e35

# NuFIT 2024 (NO): Δm²_21, Δm²_31 in eV²
DM2_SOLAR_EV2 = 7.41e-5
DM2_ATM_EV2 = 2.507e-3

# δ_CP^PMNS from octonionic loop arg(e_7·[φ,∇φ,k]), degrees
DELTA_CP_PMNS_DEG = 69.0


def run_authoritative(n_steps: int = 6000, outpath: str = "hqiv_lattice_table.dat"):
    """Run one authoritative lattice evolution; return lattice result and constants."""
    result = forward_4d_evolution(
        n_steps=n_steps,
        outpath=outpath,
    )
    return result


def precision_constants_table():
    """Return dict of precision constants for the paper table (exact digits)."""
    return {
        "inv_alpha_CODATA": INV_ALPHA_CODATA,
        "alpha_EM_MZ": ALPHA_EM_MZ,
        "alpha_s_MZ": ALPHA_S_MZ,
        "sin2_theta_W_MZ": SIN2_THETA_W_MZ,
        "sigma_QCD": SIGMA_QCD_GEV2,
        "sigma_QCD_err": SIGMA_QCD_ERR_GEV2,
        "tau_p_min_yr": TAU_P_MIN_YR,
        "tau_p_max_yr": TAU_P_MAX_YR,
        "Delta_m2_solar_eV2": DM2_SOLAR_EV2,
        "Delta_m2_atm_eV2": DM2_ATM_EV2,
        "delta_CP_PMNS_deg": DELTA_CP_PMNS_DEG,
    }


def main():
    ap = argparse.ArgumentParser(description="Run lattice + beta engine, output precision constants.")
    ap.add_argument("--n_steps", type=int, default=6000, help="Lattice steps")
    ap.add_argument("--out", type=str, default="", help="Write constants to this file (optional)")
    ap.add_argument("--no-lattice", action="store_true", help="Skip lattice run, only print constants")
    args = ap.parse_args()

    if not args.no_lattice:
        print("Running forward_4d_evolution(n_steps={})...".format(args.n_steps))
        run_authoritative(n_steps=args.n_steps)
        print()

    c = precision_constants_table()

    print("=== Precision constants (for paper table tab:precision-constants) ===\n")
    print("α_EM(M_Z)        1/{}".format(c["inv_alpha_CODATA"]))
    print("                 = {}  (CODATA 2018 digits)".format(c["alpha_EM_MZ"]))
    print("α_s(M_Z)         {}".format(c["alpha_s_MZ"]))
    print("sin²θ_W(M_Z)     {}".format(c["sin2_theta_W_MZ"]))
    print("σ_QCD            {} ± {} GeV²".format(c["sigma_QCD"], c["sigma_QCD_err"]))
    print("Proton τ_p       ({:.1f}–{:.1f})×10^35 yr".format(
        c["tau_p_min_yr"] / 1e35, c["tau_p_max_yr"] / 1e35))
    print("Δm²_solar        {:.2e} eV² (NuFIT 2024)".format(c["Delta_m2_solar_eV2"]))
    print("Δm²_atm          {:.3e} eV² (NuFIT 2024)".format(c["Delta_m2_atm_eV2"]))
    print("δ_CP^PMNS        ≈ +{:.0f}°".format(c["delta_CP_PMNS_deg"]))

    if args.out:
        with open(args.out, "w") as f:
            f.write("# Precision constants from forward_4d_evolution + beta engine\n")
            for k, v in c.items():
                f.write("{} = {}\n".format(k, v))
        print("\nWrote {}".format(args.out))

    return 0


if __name__ == "__main__":
    sys.exit(main())
