# HQIV Cosmology for CLASS / HiCLASS
# Run after installing CLASS with modified gravity support
# HiCLASS: github.com/lesgourg/class_public + HiCLASS fork for modified gravity

import numpy as np

# HQIV parameters (locked from our sandbox)
beta = 0.81          # horizon strength (a_min = beta * c * H)
omega_b = 0.0224     # visible baryons only
omega_cdm = 0.0      # ZERO cold dark matter
h = 0.70
Omega_m_vis = 0.048

# Base CLASS input dict for HQIV
# Background override: custom Hubble function
# H(a) = H0 * sqrt( Omega_r a^{-4} + Omega_m_vis * a^{-3} * (1 + delta_QI) + Omega_Lambda_QI(a) )
# where delta_QI and Omega_Lambda_QI come from our A_eff = A_std + beta * H

def get_hqiv_params():
    """Return CLASS/HiCLASS parameter dict for HQIV."""
    return {
        "output": "tCl pCl lCl mPk",
        "l_max_scalars": 3000,
        "lensing": "yes",
        "P_k_max_1/Mpc": 10,
        "omega_b": omega_b,
        "omega_cdm": 0.0,
        "H0": h * 100,
        "Omega_k": 0.0,
        "A_s": 2.1e-9,
        "n_s": 0.96,
        "tau_reio": 0.054,
        # HiCLASS modified gravity flags
        "mg_flag": 1,
        "mu0": 1.0,    # placeholder; override via custom module
        "gamma0": 1.0,
    }


def run_with_class_background_table(hqiv_Ha_path="hqiv_Ha.txt"):
    """
    Step-by-step run instructions:
    1. Run HQIV sandbox: python horizon_modes/python/bulk.py
       -> produces hqiv_Ha.txt (columns: a, H_over_H0)
    2. In CLASS source (background.c or use external table), interpolate H(a) from this table.
    3. Add to perturbation module: delta_m'' = ... + QI_boost_term (from our growth ODE).
    4. Compile CLASS/HiCLASS and run with custom background.

    This function requires CLASS Python interface (classy) and a CLASS build
    that accepts a custom H(a) table. If your CLASS build supports it:
    """
    try:
        from classy import Class
    except ImportError:
        print("CLASS Python interface not installed. Install CLASS and classy.")
        print("See class_integration/MODIFICATIONS.md for source changes.")
        return None

    cosmo = Class()
    params = get_hqiv_params()
    # Standard run (no custom background) for testing:
    cosmo.set(params)
    # If your CLASS has background_table support:
    # cosmo.set({"background_table": hqiv_Ha_path})
    cosmo.compute()
    cls = cosmo.lensed_cl(3000)
    cosmo.struct_cleanup()
    cosmo.empty()
    return cls


# Example: generate H(a) table from sandbox and save where CLASS can find it
def export_Ha_for_class(a_arr, H_arr, H0, path="hqiv_Ha.txt"):
    """Write a, H/H0 for CLASS background table."""
    H_over_H0 = np.asarray(H_arr) / H0
    np.savetxt(
        path,
        np.column_stack([np.asarray(a_arr), H_over_H0]),
        header="a  H_over_H0",
        comments="",
        fmt="%.6e",
    )


if __name__ == "__main__":
    # Generate table from sandbox (run from repo root so sandbox is on path)
    import sys
    import os
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(repo, "sandbox"))
    try:
        from hqiv_background import run_background
        from hqiv_background import H0
        a_arr, H_arr = run_background()
        export_Ha_for_class(a_arr, H_arr, H0, path=os.path.join(repo, "class_integration", "hqiv_Ha.txt"))
        print("Wrote class_integration/hqiv_Ha.txt")
    except Exception as e:
        print("Run from repo root with sandbox deps installed:", e)
