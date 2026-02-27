# CLASS / HiCLASS source modifications for HQIV

Use CLASS with modified gravity support (e.g. HiCLASS branch: [lesgourg/class_public](https://github.com/lesgourg/class_public) + HiCLASS fork).

**Authoritative HQIV→CLASS interface (paper-aligned)**

- All cosmological inputs to CLASS (background H(a), curvature, effective G(a), and fiducial parameters) must ultimately come from a **single upstream source**: the HQIV bulk/lattice generator `horizon_modes/python/bulk.py`.
- The only permitted way to change the CLASS background for HQIV is to **replace Friedmann integration by table interpolation of H(a)** built from the bulk output; do **not** introduce extra free CLASS parameters to “tune” HQIV.
- Any perturbation or low-ℓ modifications (Sections 2–3 below) must use the **same bulk-derived background** and remain consistent with the equations and choices in `paper/main.tex`. Treat them as implementation details of the HQIV framework, not as extra knobs.

## 1. Background: custom H(a)

**File:** `background.c` (or equivalent in your CLASS tree)

- **Replace** the standard Friedmann integration with **table interpolation** of the HQIV H(a) from the main bulk/lattice pipeline.
- Generate the table by running the horizon-modes code (from the HQIV repo root):
  - `python horizon_modes/python/bulk.py`  
    which writes a background / lattice table (see the top-level `README.md` and `class_hqiv_patches/README.md` for exact filenames and formats used in the paper).
- In CLASS:
  - Either add an input option `background_table = "hqiv_Ha.txt"` (or whatever name you give the extracted H(a) table) and, when set, in the background module compute H(a) by interpolating this table (and derive conformal time, etc., from it).
  - Or hardcode a wrapper that reads the table once and uses it in place of the standard H(a) from Friedmann.

## 2. Perturbations: QI growth term

**File:** `perturbations.c` (or the module that evolves δ_m)

- **Add** to the δ_m evolution (second-order equation for the matter overdensity) the QI inertia-reduction factor:
  - **Factor:** `1 / (1 − beta * c * H / |A|)` so that effective growth is boosted when the horizon term is significant.
- Implement this as a time-dependent coefficient in the growth ODE (e.g. in the source that computes δ_m'').

## 3. Horizon cutoff on long modes

- **Add** exponential damping for the longest modes:
  - `exp(−k * Θ_local / (2π))` for `k < 2π/Θ`.
- Apply this in the transfer function or in the CMB source (line-of-sight integration) for modes with k below the horizon scale so that the low-ℓ CMB power is suppressed as in the paper (≈18% at ℓ < 30).

## 4. Parameter summary for CLASS

| Parameter   | HQIV value | Note                                     |
|------------|------------|------------------------------------------|
| `omega_b`  | 0.0224     | Visible baryons only                     |
| `omega_cdm`| 0.0        | Zero dark matter                         |
| `H0`       | 70         | h = 0.70 (example; see paper fiducial)   |
| `mg_flag`  | 1          | Modified gravity on                      |
| `beta`     | —          | Deprecated here; use bulk/CLASS pipeline |

For production work aligned with the paper, **do not** use the legacy `sandbox/hqiv_background.py` ODE. Instead, always take H(a) and related background quantities from the `horizon_modes/python/bulk.py` → lattice-table → CLASS pipeline.
