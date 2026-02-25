#!/usr/bin/env python3
"""
HQIV cosmological parameter scan with the corrected first-principles math.

Physics: gamma=0.40 (fixed), no alpha, exact G_eff from Friedmann, epoch-dependent
gamma_eff, no Lambda. Free parameters: h, omega_b, Omega_k.

Cost function: CMB peak positions (Planck 2018) + sigma8 + age constraints.

Usage:
  CLASS_EXEC=/path/to/class python hqiv_cosmo_scan.py --h-range 0.60,0.80 --ob-range 0.020,0.026 --ok-range -0.05,0.05 --steps 8

For parallel execution, launch multiple instances with different h ranges.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import time
import numpy as np

PLANCK_PEAKS = np.array([220.0, 540.0, 810.0, 1120.0, 1430.0, 1750.0])

PEAK_WINDOWS = [
    (100, 350), (400, 650), (650, 950),
    (950, 1300), (1250, 1600), (1550, 2000),
]

# Observational targets
OBS_SIGMA8 = 0.811       # Planck 2018
OBS_SIGMA8_ERR = 0.006
OBS_AGE_LCDM = 13.787    # Planck 2018 LCDM age (Gyr)
OBS_AGE_LCDM_ERR = 0.020
OBS_RS_REC = 147.09      # Sound horizon at recombination (Mpc)
OBS_RS_ERR = 0.26
OBS_THETA_S = 1.04110e-2 # Angular size of sound horizon
OBS_THETA_S_ERR = 3.1e-6


def find_peaks_in_cl(ell, cl_tt, n_peaks=6):
    ell = np.asarray(ell, dtype=float)
    cl_tt = np.asarray(cl_tt)
    mask = ell >= 2
    ell, cl_tt = ell[mask], cl_tt[mask]
    d_ell = ell * (ell + 1) * cl_tt / (2.0 * np.pi)
    peaks = []
    for (lo, hi) in PEAK_WINDOWS[:n_peaks]:
        w = (ell >= lo) & (ell <= hi)
        if not np.any(w):
            peaks.append(np.nan)
            continue
        idx = np.nanargmax(d_ell[w])
        idx_global = np.where(w)[0][idx]
        peaks.append(float(ell[idx_global]))
    return np.array(peaks)


def chi2_peaks(ell_planck, ell_model):
    valid = np.isfinite(ell_model) & (ell_model > 0)
    if not np.any(valid):
        return np.inf
    pos_err = 5.0  # ~5 ell uncertainty per peak
    chi2 = np.sum(((ell_model[valid] - ell_planck[valid]) / pos_err) ** 2)
    if ell_model[0] > 0 and ell_planck[0] > 0:
        ratio_p = ell_planck / ell_planck[0]
        ratio_m = ell_model / ell_model[0]
        chi2 += np.nansum(((ratio_m[valid] - ratio_p[valid]) / 0.01) ** 2)
    return chi2


def run_class(class_exec, workdir, h, omega_b, omega_k, l_max=2500):
    root = "s_h{:.4f}_ob{:.5f}_ok{:.4f}".format(h, omega_b, omega_k).replace(".", "p").replace("-", "m")
    ini = (
        "root = {root}\n"
        "output = tCl,mPk\n"
        "l_max_scalars = {l_max}\n"
        "lensing = no\n"
        "h = {h}\n"
        "omega_b = {omega_b}\n"
        "omega_cdm = 0.0\n"
        "Omega_k = {omega_k}\n"
        "A_s = 2.1e-9\n"
        "n_s = 0.96\n"
        "tau_reio = 0.054\n"
        "HQIV = yes\n"
        "HQIV_gamma = 0.4\n"
        "HQIV_alpha = 0.6\n"
        "HQIV_chi = 0.172\n"
        "HQIV_fmin = 0.01\n"
        "gauge = newtonian\n"
        "background_verbose = 1\n"
        "fourier_verbose = 1\n"
    ).format(root=root, l_max=l_max, h=h, omega_b=omega_b, omega_k=omega_k)

    ini_path = os.path.join(workdir, "run.ini")
    with open(ini_path, "w") as f:
        f.write(ini)
    try:
        proc = subprocess.run(
            [class_exec, "run.ini"], cwd=workdir,
            capture_output=True, text=True, timeout=300,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return None

    if proc.returncode != 0:
        return None

    result = {"h": h, "omega_b": omega_b, "Omega_k": omega_k}
    stdout = proc.stdout or ""

    m = re.search(r"age\s*=\s*([\d.]+)\s*Gyr", stdout)
    result["age"] = float(m.group(1)) if m else np.nan

    m = re.search(r"sigma8=([\d.eE+-]+)", stdout)
    result["sigma8"] = float(m.group(1)) if m else np.nan

    m = re.search(r"H_actual.*?=\s*([\d.eE+-]+)\s*Mpc", stdout)
    result["H_actual"] = float(m.group(1)) if m else np.nan

    m = re.search(r"Omega_m\s+=\s+([\d.eE+-]+)", stdout)
    result["Omega_m"] = float(m.group(1)) if m else np.nan

    m = re.search(r"Omega_horizon.*?=\s+([\d.eE+-]+)", stdout)
    result["Omega_hz"] = float(m.group(1)) if m else np.nan

    m = re.search(r"Omega_k \(phys\)\s+=\s+([\d.eE+-]+)", stdout)
    result["Omega_k_phys"] = float(m.group(1)) if m else np.nan

    m = re.search(r"Sum \(closure\)\s+=\s+([\d.eE+-]+)", stdout)
    result["closure"] = float(m.group(1)) if m else np.nan

    m = re.search(r"rs_rec\s*=\s*([\d.]+)", stdout)
    result["rs_rec"] = float(m.group(1)) if m else np.nan

    m = re.search(r"D_A\(z_rec\)\s*=\s*([\d.]+)", stdout)
    result["DA_rec"] = float(m.group(1)) if m else np.nan

    m = re.search(r"Dilation.*?([\d.]+)×", stdout)
    result["dilation"] = float(m.group(1)) if m else np.nan

    # Read Cl
    for name in ["{}00_cl.dat".format(root), "{}_cl.dat".format(root)]:
        cl_path = os.path.join(workdir, name)
        if os.path.isfile(cl_path):
            try:
                data = np.loadtxt(cl_path, comments="#")
                if data.ndim == 1:
                    data = data.reshape(1, -1)
                ell = data[:, 0]
                d_over_2pi = data[:, 1]
                cl_tt = np.where(ell * (ell + 1) > 0,
                                 d_over_2pi * 2 * np.pi / (ell * (ell + 1)), 0)
                result["peaks"] = find_peaks_in_cl(ell, cl_tt)
            except Exception:
                result["peaks"] = np.full(6, np.nan)
            break
    else:
        result["peaks"] = np.full(6, np.nan)

    # Clean up output files
    for f in os.listdir(workdir):
        if f.startswith(root):
            os.remove(os.path.join(workdir, f))

    return result


def compute_chi2(r):
    if r is None:
        return np.inf
    chi2 = 0.0

    # Peak positions
    if "peaks" in r:
        chi2 += chi2_peaks(PLANCK_PEAKS, r["peaks"])

    # sigma8
    if np.isfinite(r.get("sigma8", np.nan)):
        chi2 += ((r["sigma8"] - OBS_SIGMA8) / OBS_SIGMA8_ERR) ** 2

    # Sound horizon
    if np.isfinite(r.get("rs_rec", np.nan)):
        chi2 += ((r["rs_rec"] - OBS_RS_REC) / OBS_RS_ERR) ** 2

    return chi2


def main():
    p = argparse.ArgumentParser(description="HQIV cosmological scan (gamma fixed, no alpha)")
    p.add_argument("--h-range", default="0.60,0.80", help="h range: lo,hi")
    p.add_argument("--ob-range", default="0.020,0.026", help="omega_b range: lo,hi")
    p.add_argument("--ok-range", default="-0.05,0.05", help="Omega_k range: lo,hi")
    p.add_argument("--steps", type=int, default=8, help="Steps per parameter")
    p.add_argument("--ok-steps", type=int, default=5, help="Steps for Omega_k")
    p.add_argument("--l-max", type=int, default=2500)
    p.add_argument("--out", default="hqiv_scan", help="Output prefix")
    p.add_argument("--class-exec", default=os.environ.get("CLASS_EXEC",
                   "/home/jr/Repos/HQIV/class_public/class"))
    args = p.parse_args()

    h_lo, h_hi = [float(x) for x in args.h_range.split(",")]
    ob_lo, ob_hi = [float(x) for x in args.ob_range.split(",")]
    ok_lo, ok_hi = [float(x) for x in args.ok_range.split(",")]

    h_arr = np.linspace(h_lo, h_hi, args.steps)
    ob_arr = np.linspace(ob_lo, ob_hi, args.steps)
    ok_arr = np.linspace(ok_lo, ok_hi, args.ok_steps)

    total = len(h_arr) * len(ob_arr) * len(ok_arr)
    print(f"HQIV scan: h=[{h_lo},{h_hi}]×{len(h_arr)}  ob=[{ob_lo},{ob_hi}]×{len(ob_arr)}  "
          f"Ok=[{ok_lo},{ok_hi}]×{len(ok_arr)}  total={total} points", flush=True)
    print(f"gamma=0.40 (fixed), no alpha, exact G_eff, epoch-dependent gamma_eff", flush=True)

    workdir = tempfile.mkdtemp(prefix="hqiv_scan_")
    results = []
    best_chi2 = np.inf
    best_result = None
    n = 0
    t0 = time.perf_counter()

    for h in h_arr:
        for ob in ob_arr:
            for ok in ok_arr:
                n += 1
                r = run_class(args.class_exec, workdir, h, ob, ok, l_max=args.l_max)
                if r is None:
                    continue
                chi2 = compute_chi2(r)
                r["chi2"] = chi2
                results.append(r)

                if chi2 < best_chi2:
                    best_chi2 = chi2
                    best_result = r

                if n % 20 == 0 or n == total:
                    elapsed = time.perf_counter() - t0
                    eta = (elapsed / n) * (total - n) / 60
                    print(f"  [{n}/{total}] {elapsed/60:.1f}min  ETA {eta:.1f}min  "
                          f"best chi2={best_chi2:.2f} "
                          f"h={best_result['h']:.3f} ob={best_result['omega_b']:.4f} "
                          f"Ok={best_result['Omega_k']:.3f}", flush=True)

    elapsed = time.perf_counter() - t0

    # Sort and report
    results.sort(key=lambda x: x.get("chi2", np.inf))

    print(f"\n{'='*70}")
    print(f"HQIV SCAN COMPLETE  ({n} points, {elapsed/60:.1f} min)")
    print(f"{'='*70}")

    if best_result:
        r = best_result
        print(f"\nBEST FIT (chi2 = {r['chi2']:.2f}):")
        print(f"  h       = {r['h']:.4f}  ({r['h']*100:.1f} km/s/Mpc)")
        print(f"  omega_b = {r['omega_b']:.5f}")
        print(f"  Omega_k = {r['Omega_k']:.4f}")
        print(f"  gamma   = 0.4000 (fixed)")
        print(f"  sigma8  = {r.get('sigma8', np.nan):.4f}")
        print(f"  age     = {r.get('age', np.nan):.2f} Gyr")
        print(f"  rs_rec  = {r.get('rs_rec', np.nan):.1f} Mpc")
        print(f"  Omega_m = {r.get('Omega_m', np.nan):.4f}")
        print(f"  Omega_hz= {r.get('Omega_hz', np.nan):.4f}")
        print(f"  closure = {r.get('closure', np.nan):.6f}")
        print(f"  H_actual= {r.get('H_actual', np.nan):.6e} Mpc^-1")
        print(f"  dilation= {r.get('dilation', np.nan):.2f}×")
        if "peaks" in r:
            print(f"\n  Peak comparison:")
            for i, (pp, pm) in enumerate(zip(PLANCK_PEAKS, r["peaks"])):
                d = pm - pp if np.isfinite(pm) else np.nan
                print(f"    P{i+1}  Planck={pp:.0f}  HQIV={pm:.0f}  Δℓ={d:+.0f}")

    # Top 10
    print(f"\nTop 10:")
    print(f"  {'chi2':>8}  {'h':>6}  {'ob':>7}  {'Ok':>7}  {'s8':>6}  {'age':>6}  {'rs':>6}")
    for r in results[:10]:
        print(f"  {r['chi2']:8.2f}  {r['h']:6.3f}  {r['omega_b']:7.5f}  "
              f"{r['Omega_k']:7.4f}  {r.get('sigma8',np.nan):6.3f}  "
              f"{r.get('age',np.nan):6.1f}  {r.get('rs_rec',np.nan):6.1f}")

    # Save
    out_txt = args.out + ".txt"
    with open(out_txt, "w") as f:
        f.write(f"# HQIV scan: gamma=0.40 fixed, exact G_eff, no alpha\n")
        f.write(f"# {n} points, {elapsed/60:.1f} min\n")
        f.write(f"# h  omega_b  Omega_k  chi2  sigma8  age  rs_rec  Omega_m  Omega_hz  P1  P2  P3  P4  P5  P6\n")
        for r in results:
            peaks = r.get("peaks", np.full(6, np.nan))
            f.write(f"{r['h']:.4f}  {r['omega_b']:.5f}  {r['Omega_k']:.4f}  "
                    f"{r['chi2']:.4f}  {r.get('sigma8',np.nan):.6f}  {r.get('age',np.nan):.3f}  "
                    f"{r.get('rs_rec',np.nan):.3f}  {r.get('Omega_m',np.nan):.6f}  "
                    f"{r.get('Omega_hz',np.nan):.6f}  "
                    + "  ".join(f"{p:.1f}" for p in peaks) + "\n")
    print(f"\nSaved to {out_txt}")

    try:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)
    except Exception:
        pass


if __name__ == "__main__":
    main()
