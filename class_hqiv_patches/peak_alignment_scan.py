#!/usr/bin/env python3
"""
Input-space search to best match CMB multipole peaks to Planck.

Sweeps gamma, omega_b, h, and alpha (paper: peak_alignment_scan / peak_alignment_scan3).
Uses CLASS with HQIV patches and Omega_eff closure (H(a=1)=H0).
Requires: CLASS built with class_hqiv_patches and classy Python interface.

  cd /path/to/class_public && make && export PYTHONPATH=$PWD/python:$PYTHONPATH
  python peak_alignment_scan.py [--steps 4] [--out results]

  # Fine search filling a time budget (e.g. 6 hours):
  CLASS_EXEC=/path/to/class python peak_alignment_scan.py --duration 6 --out peak_scan_6h

  # Custom fine range (gamma_lo,gamma_hi,omega_b_lo,omega_b_hi):
  CLASS_EXEC=... python peak_alignment_scan.py --duration 6 --fine-range 0.33,0.39,0.023,0.027 --out peak_scan_6h

Reference: paper/main.tex Table (Planck 2018) P1=220, P2=540, P3=810, P4=1120, P5=1430, P6=1750.
"""

from __future__ import division, print_function

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

# Planck 2018 reference positions (paper Table)
PLANCK_PEAKS = np.array([220.0, 540.0, 810.0, 1120.0, 1430.0, 1750.0])

# Search windows for finding local maxima (ell) around each Planck peak
PEAK_WINDOWS = [
    (100, 350),    # P1
    (400, 650),    # P2
    (650, 950),    # P3
    (950, 1300),   # P4
    (1250, 1600),  # P5
    (1550, 2000),  # P6
]


def find_peaks_in_cl(ell, cl_tt, n_peaks=6):
    """Return multipole positions of first n_peaks acoustic peaks (local maxima of D_ell)."""
    # D_ell = ell(ell+1)*C_ell/(2*pi)
    ell = np.asarray(ell, dtype=float)
    cl_tt = np.asarray(cl_tt)
    mask = ell >= 2
    ell = ell[mask]
    cl_tt = cl_tt[mask]
    d_ell = ell * (ell + 1) * cl_tt / (2.0 * np.pi)

    peaks_ell = []
    for (lo, hi) in PEAK_WINDOWS[:n_peaks]:
        w = (ell >= lo) & (ell <= hi)
        if not np.any(w):
            peaks_ell.append(np.nan)
            continue
        idx = np.nanargmax(d_ell[w])
        # global index
        idx_global = np.where(w)[0][idx]
        # refine: local max in window
        peaks_ell.append(float(ell[idx_global]))
    return np.array(peaks_ell)


def cost_peak_positions(ell_planck, ell_model):
    """Combined cost: position MSE + peak-ratio MSE (paper: peak-position and peak-ratio agreement)."""
    ell_planck = np.asarray(ell_planck)
    ell_model = np.asarray(ell_model)
    valid = np.isfinite(ell_model) & (ell_model > 0)
    if not np.any(valid):
        return np.inf
    # Position cost: normalized MSE
    pos_cost = np.mean((ell_model[valid] - ell_planck[valid]) ** 2) / (np.mean(ell_planck ** 2) + 1e-30)
    # Ratio cost: ell_n/ell_1 vs Planck
    if ell_model[0] > 0 and ell_planck[0] > 0:
        ratio_planck = ell_planck / ell_planck[0]
        ratio_model = ell_model / ell_model[0]
        ratio_cost = np.nanmean((ratio_model[valid] - ratio_planck[valid]) ** 2)
    else:
        ratio_cost = np.inf
    return pos_cost + ratio_cost


def _run_class_exec(class_exec, workdir, params, l_max=2500, verbose=False):
    """Run CLASS binary with generated .ini; return (ell, cl_tt, age_Gyr) or (None, None, nan)."""
    # Unique root per run so we never read another run's output (fixes identical cost across alphas)
    g = params["HQIV_gamma"]
    ob = params["omega_b"]
    h = params["h"]
    a = params["HQIV_alpha"]
    root = "scan_g{:.3f}_ob{:.4f}_h{:.3f}_a{:.3f}".format(g, ob, h, a).replace(".", "p")
    ini_content = (
        "root = {root}\n"
        "output = tCl\n"
        "l_max_scalars = {l_max}\n"
        "lensing = no\n"
        "h = {h}\n"
        "omega_b = {omega_b}\n"
        "omega_cdm = 0.0\n"
        "Omega_k = 0.0\n"
        "A_s = 2.1e-9\n"
        "n_s = 0.96\n"
        "tau_reio = 0.054\n"
        "HQIV = yes\n"
        "HQIV_gamma = {gamma}\n"
        "HQIV_alpha = {alpha}\n"
        "HQIV_chi = 0.172\n"
        "HQIV_fmin = 0.01\n"
        "gauge = newtonian\n"
        "background_verbose = 1\n"
    ).format(
        root=root,
        l_max=l_max,
        h=params["h"],
        omega_b=params["omega_b"],
        gamma=params["HQIV_gamma"],
        alpha=params["HQIV_alpha"],
    )
    ini_path = os.path.join(workdir, "run.ini")
    with open(ini_path, "w") as f:
        f.write(ini_content)
    try:
        proc = subprocess.run(
            [class_exec, "run.ini"],
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        if verbose:
            print("CLASS exec failed:", e, file=sys.stderr)
        return None, None, np.nan
    if proc.returncode != 0:
        if verbose:
            print("CLASS stderr:", proc.stderr[:500] if proc.stderr else "", file=sys.stderr)
        # Do not read existing cl file: it would be from a previous run (stale data)
        return None, None, np.nan

    # Parse age from stdout: " -> age = X.XXX Gyr"
    age_Gyr = np.nan
    if proc.stdout:
        m = re.search(r"age\s*=\s*([\d.]+)\s*Gyr", proc.stdout)
        if m:
            age_Gyr = float(m.group(1))

    # CLASS writes {root}00_cl.dat (or {root}_cl.dat): columns ell, [l(l+1)/2pi] C_l (i.e. D_l/2pi)
    for name in ["{}00_cl.dat".format(root), "{}_cl.dat".format(root)]:
        cl_path = os.path.join(workdir, name)
        if os.path.isfile(cl_path):
            break
    else:
        if verbose:
            print("No output file", os.path.join(workdir, "{}*_cl.dat".format(root)), file=sys.stderr)
        return None, None, age_Gyr
    try:
        data = np.loadtxt(cl_path, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        ell = data[:, 0]
        if data.shape[1] < 2:
            return None, None, age_Gyr
        # Column is [l(l+1)/2pi] C_l; we need C_l for find_peaks: C_l = col * 2*pi / (l*(l+1))
        d_ell_over_2pi = data[:, 1]
        cl_tt = np.where(ell * (ell + 1) > 0, d_ell_over_2pi * (2.0 * np.pi) / (ell * (ell + 1)), 0.0)
        return ell, cl_tt, age_Gyr
    except Exception as e:
        if verbose:
            print("Parse cl file failed:", e, file=sys.stderr)
        return None, None, age_Gyr


def run_class_one(params, l_max=2500, verbose=False, class_exec=None, workdir=None):
    """Run CLASS with given dict of parameters; return (ell, cl_tt, age_Gyr) or (None, None, nan) on failure."""
    if class_exec:
        return _run_class_exec(class_exec, workdir or os.getcwd(), params, l_max=l_max, verbose=verbose)

    try:
        from classy import Class
    except ImportError:
        if verbose:
            print("classy not found. Set CLASS_EXEC to path to class binary to use executable fallback.", file=sys.stderr)
        return None, None, np.nan

    base = {
        "output": "tCl",
        "l_max_scalars": int(l_max),
        "lensing": "no",
        "omega_cdm": 0.0,
        "Omega_k": 0.0,
        "A_s": 2.1e-9,
        "n_s": 0.96,
        "tau_reio": 0.054,
        "HQIV": "yes",
        "HQIV_chi": 0.172,
        "HQIV_fmin": 0.01,
    }
    base.update(params)

    try:
        cosmo = Class()
        cosmo.set(base)
        cosmo.compute()
        cl_dict = cosmo.raw_cl()
        ell = cl_dict["ell"]
        cl_tt = cl_dict["tt"]
        age_Gyr = np.nan
        if hasattr(cosmo, "get_current_derived_parameters"):
            try:
                d = cosmo.get_current_derived_parameters([])
                if isinstance(d, dict) and "age" in d:
                    age_Gyr = float(d["age"])
            except Exception:
                pass
        cosmo.struct_cleanup()
        cosmo.empty()
        return ell, cl_tt, age_Gyr
    except Exception as e:
        if verbose:
            print("CLASS run failed:", e, file=sys.stderr)
        return None, None, np.nan


def time_one_run(class_exec, workdir, l_max=2500, n_warmup=3):
    """Return mean seconds per CLASS run (n_warmup runs)."""
    params = {"HQIV_gamma": 0.36, "HQIV_alpha": 0.6, "omega_b": 0.025, "h": 0.732}
    times = []
    for _ in range(n_warmup):
        t0 = time.perf_counter()
        run_class_one(params, l_max=l_max, verbose=False, class_exec=class_exec, workdir=workdir)
        times.append(time.perf_counter() - t0)
    return float(np.mean(times))


def run_search(gamma_arr, omega_b_arr, h_arr, alpha_arr, verbose=True, class_exec=None, workdir=None, l_max=2500, progress_every=0):
    """Grid over gamma, omega_b, h, alpha; return list of (params, peaks_ell, cost, age_Gyr)."""
    results = []
    total = len(gamma_arr) * len(omega_b_arr) * len(h_arr) * len(alpha_arr)
    n = 0
    start_t = time.perf_counter()
    for gamma in gamma_arr:
        for omega_b in omega_b_arr:
            for h in h_arr:
                for alpha in alpha_arr:
                    n += 1
                    if verbose:
                        print("Run {} / {}: gamma={:.3f} omega_b={:.4f} h={:.3f} alpha={:.2f}".format(
                            n, total, gamma, omega_b, h, alpha), flush=True)
                    elif progress_every and n % progress_every == 0:
                        elapsed = time.perf_counter() - start_t
                        eta = (elapsed / n) * (total - n) if n > 0 else 0
                        print("Progress: {} / {} ({:.1f}%)  elapsed {:.1f} min  ETA {:.1f} min".format(
                            n, total, 100.0 * n / total, elapsed / 60, eta / 60), flush=True)
                    params = {
                        "HQIV_gamma": gamma,
                        "HQIV_alpha": alpha,
                        "omega_b": omega_b,
                        "h": h,
                    }
                    ell, cl_tt, age_Gyr = run_class_one(
                        params, l_max=l_max, verbose=verbose, class_exec=class_exec, workdir=workdir
                    )
                    if ell is None:
                        continue
                    peaks = find_peaks_in_cl(ell, cl_tt, n_peaks=6)
                    cost = cost_peak_positions(PLANCK_PEAKS, peaks)
                    results.append((params, peaks, cost, age_Gyr))
    return results


def main():
    p = argparse.ArgumentParser(description="Peak-alignment input-space search (CLASS-HQIV)")
    p.add_argument("--steps", type=int, default=4,
                    help="Number of steps per parameter (coarse grid)")
    p.add_argument("--fine", action="store_true", help="Run fine search around best point")
    p.add_argument("--out", type=str, default="peak_scan_results",
                    help="Output base name for .txt and .npy")
    p.add_argument("--quiet", action="store_true", help="Less stdout")
    p.add_argument("--class-exec", type=str, default=os.environ.get("CLASS_EXEC"),
                    help="Path to CLASS binary (or set CLASS_EXEC). If set, run CLASS via subprocess instead of classy.")
    p.add_argument("--duration", type=float, default=None,
                    help="Target total runtime in hours (e.g. 6). Run timing benchmark then set grid to fill this budget.")
    p.add_argument("--fine-range", type=str, default=None,
                    help="Fine grid range: gamma_lo,gamma_hi,omega_b_lo,omega_b_hi (default 0.32,0.40,0.022,0.028)")
    p.add_argument("--l-max", type=int, default=2500, help="l_max_scalars for CLASS (lower = faster, less accurate).")
    p.add_argument("--progress-every", type=int, default=0,
                    help="With --quiet, print progress every N runs (default 0). Use 500 for long runs.")
    p.add_argument("--alpha", type=float, default=None,
                    help="Fix HQIV_alpha to this value (e.g. 0.6). Default 0.60.")
    p.add_argument("--h", type=float, default=None,
                    help="Fix dimensionless Hubble h = H0/100 (e.g. 0.732). Default 0.732.")
    p.add_argument("--points", type=int, default=None,
                    help="Target number of grid points (2D gamma x omega_b). Overrides --duration when set.")
    args = p.parse_args()

    class_exec = args.class_exec
    workdir = None
    if class_exec:
        workdir = tempfile.mkdtemp(prefix="peak_scan_")
        print("Using CLASS executable: {} (workdir {})".format(class_exec, workdir), flush=True)

    h_arr = np.array([args.h if args.h is not None else 0.732])
    alpha_arr = np.array([args.alpha if args.alpha is not None else 0.60])
    l_max = args.l_max

    if args.points is not None and args.points > 0 and class_exec:
        # Fixed number of points (e.g. a few hundred for sanity checks)
        n_side = max(2, int(round(np.sqrt(args.points))))
        n_gamma = n_side
        n_omega_b = max(2, args.points // n_gamma)
        if args.fine_range:
            parts = [float(x.strip()) for x in args.fine_range.split(",")]
            g_lo, g_hi, ob_lo, ob_hi = parts[0], parts[1], parts[2], parts[3]
        else:
            g_lo, g_hi, ob_lo, ob_hi = 0.32, 0.40, 0.022, 0.028
        gamma_arr = np.linspace(g_lo, g_hi, n_gamma)
        omega_b_arr = np.linspace(ob_lo, ob_hi, n_omega_b)
        total_pts = len(gamma_arr) * len(omega_b_arr) * len(h_arr) * len(alpha_arr)
        print("Points target {} -> grid gamma {} x omega_b {} ({} points)".format(
            args.points, n_gamma, n_omega_b, total_pts), flush=True)
    elif args.duration is not None and args.duration > 0 and class_exec:
        # Duration-targeted fine search: time one run, then choose grid size to fill ~duration hours
        print("Timing benchmark ({} runs)...".format(3), flush=True)
        sec_per_point = time_one_run(class_exec, workdir, l_max=l_max, n_warmup=3)
        target_sec = args.duration * 3600.0
        # Reserve ~1% for overhead and reporting
        max_points = int(0.98 * target_sec / sec_per_point)
        max_points = max(4, max_points)
        # Fine range (centered on previous best 0.36, 0.025)
        if args.fine_range:
            parts = [float(x.strip()) for x in args.fine_range.split(",")]
            g_lo, g_hi, ob_lo, ob_hi = parts[0], parts[1], parts[2], parts[3]
        else:
            g_lo, g_hi, ob_lo, ob_hi = 0.32, 0.40, 0.022, 0.028
        # 2D grid: n_gamma * n_omega_b <= max_points, roughly square
        n_side = int(round(np.sqrt(max_points)))
        n_gamma = max(2, n_side)
        n_omega_b = max(2, max_points // n_gamma)
        gamma_arr = np.linspace(g_lo, g_hi, n_gamma)
        omega_b_arr = np.linspace(ob_lo, ob_hi, n_omega_b)
        total_pts = len(gamma_arr) * len(omega_b_arr) * len(h_arr) * len(alpha_arr)
        est_sec = total_pts * sec_per_point
        print("Target {:.1f} h -> ~{} points at {:.2f} s/point -> grid gamma {} x omega_b {} ({} points, est. {:.2f} h)".format(
            args.duration, max_points, sec_per_point, n_gamma, n_omega_b, total_pts, est_sec / 3600), flush=True)
    else:
        # Fixed steps (original behavior)
        gamma_arr = np.linspace(0.36, 0.44, max(2, args.steps))
        omega_b_arr = np.linspace(0.025, 0.030, max(2, args.steps))

    print("Running peak-alignment search (CLASS-HQIV with Omega_eff closure)...", flush=True)
    print("Planck reference peaks (ell):", PLANCK_PEAKS.tolist(), flush=True)
    progress_every = args.progress_every if args.quiet else 0
    if args.duration and args.quiet and progress_every == 0:
        progress_every = 500  # default progress when running long quiet job
    run_start = time.perf_counter()
    results = run_search(
        gamma_arr, omega_b_arr, h_arr, alpha_arr,
        verbose=not args.quiet, class_exec=class_exec, workdir=workdir, l_max=l_max,
        progress_every=progress_every,
    )
    run_elapsed = time.perf_counter() - run_start
    if args.duration and run_elapsed > 0:
        print("Search completed in {:.1f} min ({:.2f} h)".format(run_elapsed / 60, run_elapsed / 3600), flush=True)

    if not results:
        print("No successful CLASS runs. Check classy and CLASS build.", file=sys.stderr)
        sys.exit(1)

    # Sort by cost
    results.sort(key=lambda x: x[2])
    best_params, best_peaks, best_cost, best_age = results[0]

    # Summary
    print("\n--- Best match (combined peak position + ratio cost) ---")
    print("  gamma={:.4f} omega_b={:.4f} h={:.3f} alpha={:.2f}".format(
        best_params["HQIV_gamma"], best_params["omega_b"], best_params["h"], best_params["HQIV_alpha"]))
    print("  cost = {:.6f}".format(best_cost))
    if not np.isnan(best_age):
        print("  age  = {:.2f} Gyr".format(best_age))
    print("\n  Peak comparison (Planck vs HQIV):")
    for i, (pp, ph) in enumerate(zip(PLANCK_PEAKS, best_peaks)):
        d = ph - pp if np.isfinite(ph) else np.nan
        print("    P{}  Planck={:.0f}  HQIV={:.0f}  Delta_ell={:+.0f}".format(i + 1, pp, ph, d))

    # Save
    out_txt = args.out + ".txt"
    out_npy = args.out + ".npy"
    with open(out_txt, "w") as f:
        f.write("# peak_alignment_scan best and full results\n")
        f.write("best_params gamma={} omega_b={} h={} alpha={}\n".format(
            best_params["HQIV_gamma"], best_params["omega_b"], best_params["h"], best_params["HQIV_alpha"]))
        f.write("best_cost={} best_age_Gyr={}\n".format(best_cost, best_age))
        f.write("Planck peaks: {}\n".format(PLANCK_PEAKS.tolist()))
        f.write("HQIV  peaks: {}\n".format(best_peaks.tolist()))
        if args.duration:
            f.write("duration_h={} n_points={}\n".format(args.duration, len(results)))
    np.save(out_npy, {"results": results, "planck_peaks": PLANCK_PEAKS, "best_params": best_params})
    print("\nWrote {} and {}".format(out_txt, out_npy))

    if workdir and os.path.isdir(workdir):
        try:
            import shutil
            shutil.rmtree(workdir, ignore_errors=True)
        except Exception:
            pass

    if args.fine and args.steps <= 4 and args.duration is None:
        print("\n--- Fine search around best ---")
        g0, ob0 = best_params["HQIV_gamma"], best_params["omega_b"]
        gamma_fine = np.linspace(max(0.3, g0 - 0.03), min(0.5, g0 + 0.03), 5)
        omega_b_fine = np.linspace(max(0.02, ob0 - 0.003), min(0.035, ob0 + 0.003), 5)
        fine_workdir = tempfile.mkdtemp(prefix="peak_scan_fine_") if class_exec else None
        fine_results = run_search(
            gamma_fine, omega_b_fine, h_arr, alpha_arr,
            verbose=not args.quiet, class_exec=class_exec, workdir=fine_workdir, l_max=l_max,
            progress_every=0,
        )
        if fine_workdir and os.path.isdir(fine_workdir):
            try:
                shutil.rmtree(fine_workdir, ignore_errors=True)
            except Exception:
                pass
        fine_results.sort(key=lambda x: x[2])
        if fine_results:
            fp, fpeaks, fc, fa = fine_results[0]
            print("  Fine best: gamma={:.4f} omega_b={:.4f} cost={:.6f} age={:.2f} Gyr".format(
                fp["HQIV_gamma"], fp["omega_b"], fc, fa))
            np.save(args.out + "_fine.npy", {"results": fine_results})

    return 0


if __name__ == "__main__":
    sys.exit(main())
