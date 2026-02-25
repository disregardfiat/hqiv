#!/usr/bin/env python3
"""
Run 200 CLASS-HQIV runs over Omega_m = 0.02 to 0.04 (step 0.0001) with T_cmb = 2.725.
Cost = multipole delta to current CMB observations (Planck 2018 peak positions).

In HQIV baryon-only, Omega_m = Omega_b, so omega_b = Omega_m * h^2.

Usage:
  CLASS_EXEC=/path/to/class python3 background_cost_scan.py [--dry-run] [--out results.csv]
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

# Planck 2018 reference peak positions (paper / peak_alignment_scan)
PLANCK_PEAKS = np.array([220.0, 540.0, 810.0, 1120.0, 1430.0, 1750.0])

PEAK_WINDOWS = [
    (100, 350), (400, 650), (650, 950),
    (950, 1300), (1250, 1600), (1550, 2000),
]


def find_peaks_in_cl(ell, cl_tt, n_peaks=6):
    """Return multipole positions of first n_peaks acoustic peaks (local maxima of D_ell)."""
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


def cost_peak_positions(ell_planck, ell_model):
    """Cost = multipole delta to observations: position MSE + peak-ratio MSE."""
    ell_planck = np.asarray(ell_planck)
    ell_model = np.asarray(ell_model)
    valid = np.isfinite(ell_model) & (ell_model > 0)
    if not np.any(valid):
        return np.inf
    pos_cost = np.mean((ell_model[valid] - ell_planck[valid]) ** 2) / (np.mean(ell_planck ** 2) + 1e-30)
    if ell_model[0] > 0 and ell_planck[0] > 0:
        ratio_planck = ell_planck / ell_planck[0]
        ratio_model = ell_model / ell_model[0]
        ratio_cost = np.nanmean((ratio_model[valid] - ratio_planck[valid]) ** 2)
    else:
        ratio_cost = np.inf
    return pos_cost + ratio_cost


def run_one(class_exec, workdir, Omega_m, h=0.73, T_cmb=2.725, l_max=2500):
    """Run CLASS once; return Omega_m, ok, cost (multipole delta), age_Gyr, peaks_ell, time_s."""
    omega_b = Omega_m * (h ** 2)
    root = "om{:05d}".format(int(round(Omega_m * 100000)))
    ini = (
        "root = {root}\n"
        "output = tCl\n"
        "l_max_scalars = {l_max}\n"
        "lensing = no\n"
        "h = {h}\n"
        "omega_b = {omega_b}\n"
        "omega_cdm = 0.0\n"
        "Omega_k = 0.0\n"
        "T_cmb = {T_cmb}\n"
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
        "thermodynamics_verbose = 0\n"
        "perturbations_verbose = 0\n"
        "fourier_verbose = 0\n"
    ).format(root=root, h=h, omega_b=omega_b, T_cmb=T_cmb, l_max=l_max)
    ini_path = os.path.join(workdir, "run.ini")
    with open(ini_path, "w") as f:
        f.write(ini)
    try:
        t0 = time.perf_counter()
        proc = subprocess.run(
            [class_exec, "run.ini"],
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=300,
        )
        elapsed = time.perf_counter() - t0
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        return {"Omega_m": Omega_m, "ok": False, "cost": np.nan, "age_Gyr": np.nan, "peaks": None, "time_s": None, "error": str(e)}

    if proc.returncode != 0:
        return {"Omega_m": Omega_m, "ok": False, "cost": np.nan, "age_Gyr": np.nan, "peaks": None, "time_s": elapsed, "stderr": (proc.stderr or "")[:500]}

    age_Gyr = np.nan
    if proc.stdout:
        # CLASS prints " -> age = X.XX Gyr" (requires background_verbose >= 1)
        m = re.search(r"age\s*=\s*([\d.]+)\s*Gyr", proc.stdout, re.IGNORECASE)
        if m:
            age_Gyr = float(m.group(1))

    # Parse Cl file: CLASS writes root00_cl.dat (or root_cl.dat)
    for name in ["{}00_cl.dat".format(root), "{}_cl.dat".format(root)]:
        cl_path = os.path.join(workdir, name)
        if os.path.isfile(cl_path):
            break
    else:
        return {"Omega_m": Omega_m, "ok": False, "cost": np.nan, "age_Gyr": age_Gyr, "peaks": None, "time_s": elapsed, "error": "No cl file"}

    try:
        data = np.loadtxt(cl_path, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        ell = data[:, 0]
        if data.shape[1] < 2:
            return {"Omega_m": Omega_m, "ok": False, "cost": np.nan, "age_Gyr": age_Gyr, "peaks": None, "time_s": elapsed}
        d_ell_over_2pi = data[:, 1]
        cl_tt = np.where(ell * (ell + 1) > 0, d_ell_over_2pi * (2.0 * np.pi) / (ell * (ell + 1)), 0.0)
    except Exception as e:
        return {"Omega_m": Omega_m, "ok": False, "cost": np.nan, "age_Gyr": age_Gyr, "peaks": None, "time_s": elapsed, "error": str(e)}

    peaks = find_peaks_in_cl(ell, cl_tt, n_peaks=6)
    cost = cost_peak_positions(PLANCK_PEAKS, peaks)
    return {"Omega_m": Omega_m, "ok": True, "cost": cost, "age_Gyr": age_Gyr, "peaks": peaks, "time_s": elapsed}


def main():
    p = argparse.ArgumentParser(description="Omega_m 0.02–0.04 scan; cost = multipole delta to Planck")
    p.add_argument("--class-exec", type=str, default=os.environ.get("CLASS_EXEC"), help="Path to CLASS binary")
    p.add_argument("--Omega-m-min", type=float, default=0.02, help="Minimum Omega_m")
    p.add_argument("--Omega-m-max", type=float, default=0.04, help="Maximum Omega_m")
    p.add_argument("--step", type=float, default=0.0001, help="Step (default 0.0001 => 200 points)")
    p.add_argument("--h", type=float, default=0.73, help="Fixed h for omega_b = Omega_m * h^2")
    p.add_argument("--T_cmb", type=float, default=2.725, help="Target T_cmb (K)")
    p.add_argument("--l-max", type=int, default=2500, help="l_max_scalars for C_l")
    p.add_argument("--dry-run", action="store_true", help="Only print run count and exit")
    p.add_argument("--out", type=str, default=None, help="Write Omega_m, cost, age, peaks to CSV")
    args = p.parse_args()

    class_exec = args.class_exec
    if not class_exec or not os.path.isfile(class_exec):
        print("Set CLASS_EXEC to path to class binary.", file=sys.stderr)
        sys.exit(1)

    Omega_m_vals = []
    x = args.Omega_m_min
    while x <= args.Omega_m_max + 0.5 * args.step:
        Omega_m_vals.append(round(x, 6))
        x += args.step
    n_runs = len(Omega_m_vals)

    if args.dry_run:
        print("Dry run: {} CLASS jobs (Omega_m {} to {} step {}), cost = multipole delta to Planck".format(
            n_runs, args.Omega_m_min, args.Omega_m_max, args.step))
        return

    workdir = tempfile.mkdtemp(prefix="class_bg_scan_")
    print("Workdir: {}".format(workdir), flush=True)
    print("Running {} runs (Omega_m [{}, {}], step {}), T_cmb = {} K, cost = multipole delta to Planck ...".format(
        n_runs, args.Omega_m_min, args.Omega_m_max, args.step, args.T_cmb), flush=True)

    total_start = time.perf_counter()
    results = []
    for i, Om in enumerate(Omega_m_vals):
        r = run_one(class_exec, workdir, Om, h=args.h, T_cmb=args.T_cmb, l_max=args.l_max)
        results.append(r)
        if not r["ok"]:
            print("  [{}] Omega_m={:.4f} FAILED: {}".format(i + 1, Om, r.get("error", r.get("stderr", "?"))[:60]), flush=True)
        elif (i + 1) % 50 == 0 or i == 0:
            print("  [{}]/{} Omega_m={:.4f} cost={:.6f}".format(i + 1, n_runs, Om, r["cost"]), flush=True)
    total_elapsed = time.perf_counter() - total_start

    ok = [r for r in results if r["ok"]]
    costs = [r["cost"] for r in ok]
    if not costs:
        print("No successful runs.", file=sys.stderr)
        sys.exit(1)

    best = min(ok, key=lambda r: r["cost"])
    print("")
    print("--- Cost (multipole delta to current observations) ---")
    print("  Successful runs: {} / {}".format(len(ok), n_runs))
    print("  Best Omega_m     = {:.4f}  (cost = {:.6f})".format(best["Omega_m"], best["cost"]))
    if best.get("peaks") is not None:
        print("  Model peaks (ell)= {}".format(" ".join("{:.0f}".format(p) for p in best["peaks"])))
        print("  Planck peaks     = {}".format(" ".join("{:.0f}".format(p) for p in PLANCK_PEAKS)))
    print("  Age at best      = {:.2f} Gyr".format(best["age_Gyr"]))
    print("  Cost range       = [{:.6f}, {:.6f}]".format(min(costs), max(costs)))
    print("  Wall time        = {:.1f} s".format(total_elapsed))

    if args.out and results:
        with open(args.out, "w") as f:
            f.write("Omega_m,cost,age_Gyr,ok,P1,P2,P3,P4,P5,P6\n")
            for r in results:
                peaks = r.get("peaks")
                p1p6 = ",".join("{:.1f}".format(peaks[i]) if peaks is not None and i < len(peaks) and np.isfinite(peaks[i]) else "" for i in range(6))
                f.write("{:.6f},{},{},{},{}\n".format(
                    r["Omega_m"],
                    r["cost"] if np.isfinite(r["cost"]) else "",
                    r["age_Gyr"] if np.isfinite(r["age_Gyr"]) else "",
                    1 if r["ok"] else 0,
                    p1p6))
        print("  Wrote {}".format(args.out))

    try:
        for f in os.listdir(workdir):
            os.unlink(os.path.join(workdir, f))
        os.rmdir(workdir)
    except OSError:
        pass


if __name__ == "__main__":
    main()
