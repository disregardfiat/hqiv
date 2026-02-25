#!/usr/bin/env python3
"""
Run a single CLASS-HQIV run with full output (tCl + mPk) and report sigma8 and
all observables that CLASS prints, for comparison to observations.

Usage:
  CLASS_EXEC=/path/to/class python3 run_single_full_output.py [--gamma 0.4] [--omega_b 0.018] [--h 0.734] [--alpha 0.6]
  Output: stdout summary + optional --out report.txt (and full CLASS log to report.log).
"""

from __future__ import division, print_function

import argparse
import os
import re
import subprocess
import sys
import tempfile

def main():
    p = argparse.ArgumentParser(description="Single CLASS-HQIV run with full output (sigma8, etc.)")
    p.add_argument("--gamma", type=float, default=0.4, help="HQIV_gamma")
    p.add_argument("--omega_b", type=float, default=0.018, help="omega_b")
    p.add_argument("--h", type=float, default=0.734, help="h = H0/100")
    p.add_argument("--alpha", type=float, default=0.6, help="HQIV_alpha")
    p.add_argument("--class-exec", type=str, default=os.environ.get("CLASS_EXEC"), help="Path to CLASS binary")
    p.add_argument("--out", type=str, default=None, help="Base name for report .txt and .log (default: print only)")
    p.add_argument("--l-max", type=int, default=2500, help="l_max_scalars")
    args = p.parse_args()

    class_exec = args.class_exec
    if not class_exec or not os.path.isfile(class_exec):
        print("Set CLASS_EXEC to path to class binary.", file=sys.stderr)
        sys.exit(1)

    root = "valley_run"
    ini_content = (
        "root = {root}\n"
        "output = tCl,mPk\n"
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
        "background_verbose = 2\n"
        "thermodynamics_verbose = 1\n"
        "fourier_verbose = 1\n"
    ).format(
        root=root,
        l_max=args.l_max,
        h=args.h,
        omega_b=args.omega_b,
        gamma=args.gamma,
        alpha=args.alpha,
    )

    workdir = tempfile.mkdtemp(prefix="class_full_")
    ini_path = os.path.join(workdir, "run.ini")
    with open(ini_path, "w") as f:
        f.write(ini_content)

    try:
        proc = subprocess.run(
            [class_exec, "run.ini"],
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        print("CLASS run failed:", e, file=sys.stderr)
        sys.exit(1)

    stdout_text = proc.stdout or ""
    stderr_text = proc.stderr or ""

    if proc.returncode != 0:
        print("CLASS exited with code", proc.returncode, file=sys.stderr)
        print("stderr:", stderr_text[:2000], file=sys.stderr)
        sys.exit(1)

    # Parse key quantities from stdout (CLASS prints " -> key = value" or " -> key=value")
    def find_float(pat, text, default=None):
        m = re.search(pat, text)
        return float(m.group(1)) if m else default

    def find_line(pat, text):
        m = re.search(pat, text, re.MULTILINE)
        return m.group(1).strip() if m else None

    age_Gyr = find_float(r"->\s*age\s*=\s*([\d.]+)\s*Gyr", stdout_text)
    conformal_age_Mpc = find_float(r"->\s*conformal age\s*=\s*([\d.]+)\s*Mpc", stdout_text)
    N_eff = find_float(r"->\s*N_eff\s*=\s*([\d.e+-]+)", stdout_text)
    z_eq = find_float(r"radiation/matter equality at z\s*=\s*([\d.]+)", stdout_text)
    sigma8_m = find_float(r"->\s*sigma8=([\d.eE+-]+)\s+for total matter", stdout_text)
    sigma8_cb = find_float(r"sigma8=([\d.eE+-]+)\s+for baryons\+cdm", stdout_text)
    z_rec = find_float(r"recombination.*at z\s*=\s*([\d.]+)", stdout_text)
    rs_rec = find_float(r"comoving sound horizon\s*=\s*([\d.]+)\s*Mpc", stdout_text)
    theta_s = find_float(r"sound horizon angle 100\*theta_s\s*=\s*([\d.]+)", stdout_text)
    theta_star = find_float(r"angle 100\*theta_\*\s*=\s*([\d.]+)", stdout_text)
    z_reio = find_float(r"reionization at z\s*=\s*([\d.]+)", stdout_text)
    # Omega from budget (look for Cosmological Constant line; CLASS may print "Bayrons" typo)
    Omega_Lambda = find_float(r"Cosmological Constant\s+Omega\s*=\s*([\d.e+-]+)", stdout_text)
    Omega_b = find_float(r"Bayrons?\s+Omega\s*=\s*([\d.e+-]+)", stdout_text)
    if Omega_b is None:
        Omega_b = find_float(r"Baryons\s+Omega\s*=\s*([\d.e+-]+)", stdout_text)
    # theta_s: "sound horizon angle 100*theta_s = X"
    if theta_s is None:
        theta_s = find_float(r"100\*theta_s\s*=\s*([\d.]+)", stdout_text)

    # S8 = sigma8 * sqrt(Omega_m/0.3); in HQIV baryons-only Omega_m ≈ Omega_b
    Omega_m = Omega_b if Omega_b is not None else (args.omega_b / (args.h**2))
    S8 = sigma8_m * (Omega_m / 0.3) ** 0.5 if sigma8_m is not None and Omega_m else None

    # Build report
    lines = []
    lines.append("=" * 60)
    lines.append("CLASS-HQIV single run: full observables (valley point)")
    lines.append("=" * 60)
    lines.append("Input: gamma={:.4f} omega_b={:.4f} h={:.4f} alpha={:.2f}".format(
        args.gamma, args.omega_b, args.h, args.alpha))
    lines.append("")
    lines.append("--- Derived (from CLASS stdout) ---")
    if age_Gyr is not None:
        lines.append("  age                    = {:.4f} Gyr".format(age_Gyr))
    if conformal_age_Mpc is not None:
        lines.append("  conformal age           = {:.4f} Mpc".format(conformal_age_Mpc))
    if N_eff is not None:
        lines.append("  N_eff                  = {:.4f}".format(N_eff))
    if z_eq is not None:
        lines.append("  z_eq                   = {:.2f}".format(z_eq))
    if Omega_Lambda is not None:
        lines.append("  Omega_Lambda (eff)      = {:.4f}".format(Omega_Lambda))
    if Omega_b is not None:
        lines.append("  Omega_b                = {:.6f}".format(Omega_b))
    lines.append("")
    if sigma8_m is not None:
        lines.append("  sigma8 (total matter)   = {:.6g}".format(sigma8_m))
        if sigma8_m > 2 or sigma8_m < 0.01:
            lines.append("    (WARNING: value far from Planck ~0.81; check HQIV P(k)/sigma8 implementation)")
    if sigma8_cb is not None:
        lines.append("  sigma8 (baryons+cdm)    = {:.6g}".format(sigma8_cb))
    if S8 is not None and 0.01 <= sigma8_m <= 2:
        lines.append("  S8 = sigma8*sqrt(Om/0.3)= {:.6f}".format(S8))
    lines.append("")
    if z_rec is not None:
        lines.append("  z_rec                  = {:.2f}".format(z_rec))
    if rs_rec is not None:
        lines.append("  sound horizon at rec    = {:.4f} Mpc".format(rs_rec))
    if theta_s is not None:
        lines.append("  100*theta_s            = {:.6f}".format(theta_s))
    if theta_star is not None:
        lines.append("  100*theta_*             = {:.6f}".format(theta_star))
    if z_reio is not None:
        lines.append("  z_reio                 = {:.2f}".format(z_reio))
    lines.append("")
    lines.append("--- Comparison to observations (reference values) ---")
    lines.append("  Planck 2018 (LambdaCDM): sigma8 ~ 0.811, S8 ~ 0.834, age ~ 13.8 Gyr")
    lines.append("  DES/Y3 weak lensing:    S8 ~ 0.76-0.78")
    lines.append("  100*theta_s (Planck):   ~ 1.041")
    lines.append("  rs (sound horizon):     ~ 147 Mpc (Planck)")
    if sigma8_m is not None and 0.01 <= sigma8_m <= 2:
        lines.append("  Delta sigma8 (HQIV - 0.81) = {:.4f}".format(sigma8_m - 0.811))
    if S8 is not None and 0.01 <= S8 <= 2:
        lines.append("  Delta S8 (HQIV - 0.834)    = {:.4f}".format(S8 - 0.834))
    if age_Gyr is not None:
        lines.append("  Delta age (HQIV - 13.8) Gyr = {:.2f}".format(age_Gyr - 13.8))
    lines.append("=" * 60)

    report = "\n".join(lines)
    print(report)

    if args.out:
        with open(args.out + ".txt", "w") as f:
            f.write(report)
        with open(args.out + ".log", "w") as f:
            f.write(stdout_text)
            if stderr_text:
                f.write("\n--- stderr ---\n")
                f.write(stderr_text)
        print("Wrote {} and {}".format(args.out + ".txt", args.out + ".log"), file=sys.stderr)

    # Optional: list output files
    for ext in ["00_cl.dat", "00_pk.dat"]:
        path = os.path.join(workdir, root + ext)
        if os.path.isfile(path):
            size = os.path.getsize(path)
            print("Output file {} ({:.1f} kB)".format(root + ext, size / 1024), file=sys.stderr)

    try:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)
    except Exception:
        pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
