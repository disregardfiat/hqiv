# CLASS / HiCLASS source modifications for HQIV

Use CLASS with modified gravity support (e.g. HiCLASS branch: [lesgourg/class_public](https://github.com/lesgourg/class_public) + HiCLASS fork).

## 1. Background: custom H(a)

**File:** `background.c` (or equivalent in your CLASS tree)

- **Replace** the standard Friedmann integration with **table interpolation** of our solved H(a) from the quadratic A_eff solver.
- Our sandbox writes `hqiv_Ha.txt` with columns: `a`, `H_over_H0`.
- In CLASS:
  - Either add an input option `background_table = "hqiv_Ha.txt"` and, when set, in the background module compute H(a) by interpolating this table (and derive conformal time, etc., from it).
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

| Parameter   | HQIV value | Note                    |
|------------|------------|-------------------------|
| `omega_b`  | 0.0224     | Visible baryons only    |
| `omega_cdm`| 0.0        | Zero dark matter        |
| `H0`       | 70         | h = 0.70                |
| `mg_flag`  | 1          | Modified gravity on     |
| `beta`     | 0.81       | Horizon strength (custom) |

Generate the background table by running from the repo root:

```bash
python sandbox/hqiv_background.py
# -> writes sandbox/hqiv_Ha.txt; copy to class_integration/ or point CLASS to it
```

Then point CLASS’s custom background to `hqiv_Ha.txt` (or the path you use).
