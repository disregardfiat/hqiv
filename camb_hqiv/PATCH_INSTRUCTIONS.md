# HQIV (Horizon-Quantized Informational Vacuum) Patch for CAMB

This directory contains patches and scripts to add **HQIV cosmology** (derived from Quantised Inertia / MiHsC) to the public [CAMB](https://github.com/cmbant/CAMB) code.

## Prerequisites

- **Fresh clone of CAMB:**  
  `git clone https://github.com/cmbant/CAMB.git && cd CAMB`
- **Fortran compiler:** gfortran 6+ (or Intel ifort)
- **Python 3.8+** and `pip` (for wrapper and tests)
- **HQIV background table:** `hqiv_Ha.txt` (columns: `a`, `H_over_H0`), e.g. from your solver in `sandbox/hqiv_Ha.txt`

## Step-by-step: apply the patch

1. **Clone CAMB and go to the Fortran tree**
   ```bash
   git clone https://github.com/cmbant/CAMB.git
   cd CAMB/fortran
   ```

2. **Apply the unified patch (Fortran type/param changes only)**
   ```bash
   patch -p1 < /path/to/HQIV/camb_hqiv/patches/camb_hqiv.patch
   ```
   If the patch fails (e.g. line offsets), apply the edits by hand using the snippets in `patches/`:
   - `hqiv_results_snippet.f90` – LoadHQIVTable, GetHoverH0, GetGratio, SetParams/Hofz/low-ℓ
   - `hqiv_equations_snippet.f90` – dtauda HQIV branch
   - `hqiv_camb_readparams_snippet.f90` – read HQIV and options in CAMB_ReadParams
   - `hqiv_python_model_snippet.py` – add HQIV fields to Python `CAMBparams._fields_` (same order as Fortran)

3. **Copy the HQIV background table**
   ```bash
   cp /path/to/HQIV/sandbox/hqiv_Ha.txt ./
   ```
   Or set the path in your `params.ini` (see below).

4. **Build CAMB**
   ```bash
   make clean
   make
   ```
   Or build the Python extension from the CAMB repo root:
   ```bash
   cd ..
   pip install -e .
   ```

5. **Run the test script**
   From the CAMB repo root (or with `PYTHONPATH` pointing to the built `camb`):
   ```bash
   python /path/to/HQIV/camb_hqiv/run_hqiv_test.py
   ```
   Or run the Python wrapper (see below) to produce the TT spectrum and plot.

## What gets modified

- **fortran/model.f90** – Add `HQIV` flag and `hqiv_Ha_file` to `CAMBparams`; add `hqiv_cs2_fac`, `hqiv_l_cut`, `hqiv_beta`.
- **fortran/results.f90** – HQIV table storage in `CAMBdata`, load `hqiv_Ha.txt`, `GetHoverH0`, `GetGratio`; in `SetParams` set `omch2=0` and flat universe when HQIV; low-ℓ damping when writing Cl.
- **fortran/equations.f90** – `dtauda` uses tabulated H(a) when HQIV; optional sound-speed factor and G(a) in perturbation sources (stubbed where needed).
- **fortran/camb.f90** – Read `HQIV`, `hqiv_Ha_file`, and optional HQIV params from the `.ini` in `CAMB_ReadParams`.
- **inifiles/** – Example `params_hqiv.ini` for HQIV runs.

When `HQIV = F` (default), behaviour is unchanged. When `HQIV = T`, the code uses the tabulated background, optional varying G(a), optional sound-speed factor, and low-ℓ cutoff as described in FILE_LIST.md and the source comments.

## params.ini options (HQIV block)

```ini
HQIV = T
hqiv_Ha_file = hqiv_Ha.txt
# Optional (defaults in code):
# hqiv_cs2_fac = 1.0
# hqiv_l_cut = 38
# hqiv_beta = 0.78
```

- **hqiv_Ha_file** – Path to table with columns `a`, `H_over_H0` (relative to run directory or absolute).
- **hqiv_cs2_fac** – Placeholder factor for effective sound speed (~0.95–1.05); 1 = standard.
- **hqiv_l_cut** – Scale for low-ℓ damping: suppression exp(-(ℓ/ℓ_cut)^1.8) for ℓ < ~80.
- **hqiv_beta** – β for horizon motive (used in G(a) and documentation; default 0.78).

## Validation

- **Compile:** `make` in `fortran/` must complete without errors.
- **Run:** `python run_hqiv_test.py` (or `camb params_hqiv.ini` from `fortran/`) runs without errors.
- **Output:** Test script prints universe age today, ℓ of first acoustic peak, and low-ℓ suppression (ℓ=2–30).

## Files in this bundle

- `PATCH_INSTRUCTIONS.md` (this file)
- `FILE_LIST.md` – List of modified files and roles
- `patches/camb_hqiv.patch` – Unified diff for all changes
- `patches/` – Optional per-file patches
- `inifiles/params_hqiv.ini` – Example HQIV params
- `run_hqiv_test.py` – Test script (age, first peak, low-ℓ suppression)
- `hqiv_camb_wrapper.py` – Short Python wrapper to run HQIV and save TT plot
