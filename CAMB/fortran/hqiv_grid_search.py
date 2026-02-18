#!/usr/bin/env python3
"""
HQIV grid search over hqiv_coupling_strength.
Run from CAMB/fortran; each CAMB run ~12 min; 12 points ≈ 2.5 h.
Results saved to hqiv_grid_results.csv. Review table in the morning.

  cd CAMB/fortran && nohup python3 hqiv_grid_search.py > hqiv_grid_stdout.log 2>&1 &
  # next morning: cat hqiv_grid_results.csv
"""
import subprocess
import re
import sys
from pathlib import Path

# Run from script directory (CAMB/fortran)
SCRIPT_DIR = Path(__file__).resolve().parent


def parse_log(log_path):
    with open(log_path, "r") as f:
        log = f.read()
    age = re.search(r"Universe age today\s*=\s*([\d.]+)", log)
    sigma8 = re.search(r"sigma_8 \(z=0\)\s*=\s*([\d.]+)", log)
    peak = re.search(r"First acoustic peak ell\s*=\s*([\d.]+)", log)
    # Log format: "  Low-ell suppression      =  -1.7984E+03 %" or "  2.0500E+01 %"
    supp = re.search(r"Low-ell suppression\s*=\s*([-\d.E+]+)", log)
    return {
        "age": float(age.group(1)) if age else None,
        "sigma8": float(sigma8.group(1)) if sigma8 else None,
        "peak": float(peak.group(1)) if peak else None,
        "supp": float(supp.group(1)) if supp else None,
    }


def main():
    try:
        import numpy as np
    except ImportError:
        np = None

    # Grid of coupling strengths: 12 runs ≈ 2.5 hours
    strengths = [round(x, 2) for x in [1.5 + (5.0 - 1.5) * i / 11 for i in range(12)]]
    results = []

    print("HQIV grid search starting... (will run ~2.5 h for 12 points)")
    print(f"Working directory: {SCRIPT_DIR}")
    print(f"Strengths: {strengths}")

    for val in strengths:
        print(f"\n=== Trial {val:.2f} ===")
        sys.stdout.flush()

        ini = SCRIPT_DIR / "params_hqiv_covariant_test.ini"
        lines = ini.read_text().splitlines()
        new_lines = []
        for line in lines:
            if line.strip().startswith("hqiv_coupling_strength"):
                new_lines.append(f"hqiv_coupling_strength = {val:.2f}")
            else:
                new_lines.append(line)
        ini.write_text("\n".join(new_lines) + "\n")

        log_path = SCRIPT_DIR / f"run_{val:.2f}.log"
        with open(log_path, "w") as logfile:
            ret = subprocess.run(
                ["stdbuf", "-oL", "./camb", "params_hqiv_covariant_test.ini"],
                cwd=SCRIPT_DIR,
                stdout=logfile,
                stderr=subprocess.STDOUT,
            )
        if ret.returncode != 0:
            print(f"  WARNING: CAMB exited with code {ret.returncode}")
        data = parse_log(log_path)
        data["strength"] = val
        results.append(data)
        print(
            f"  Age: {data['age']} Gyr | sigma8: {data['sigma8']} | Peak: {data['peak']} | Low-l supp: {data['supp']}"
        )
        sys.stdout.flush()

    csv_path = SCRIPT_DIR / "hqiv_grid_results.csv"
    with open(csv_path, "w") as f:
        f.write("strength,age,sigma8,peak,low_l_supp\n")
        for r in results:
            f.write(
                f"{r['strength']:.2f},{r['age'] or ''},{r['sigma8'] or ''},{r['peak'] or ''},{r['supp'] or ''}\n"
            )

    print("\n=== GRID SEARCH COMPLETE ===")
    print(f"Results saved to {csv_path}")
    print("Best runs will be obvious in the morning.")


if __name__ == "__main__":
    main()
