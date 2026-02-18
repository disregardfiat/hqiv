# HQIV modification – file list

All paths relative to the **CAMB repository root** (so Fortran sources live under `fortran/`).

## Files to modify

| File | Role |
|------|------|
| **fortran/model.f90** | Add to `CAMBparams`: `HQIV` (logical), `hqiv_Ha_file` (path), `hqiv_cs2_fac`, `hqiv_l_cut`, `hqiv_beta`. |
| **fortran/results.f90** | HQIV table in `CAMBdata`; load `hqiv_Ha.txt` in `SetParams`; `GetHoverH0(this,a)`, `GetGratio(this,a)`; when HQIV set `omch2=0` and flat; apply low-ℓ damping to Cl for ℓ < ~80. |
| **fortran/equations.f90** | In `dtauda`: if HQIV use tabulated H(a) (via `GetHoverH0`); in perturbation equations use optional G(a) and cs2 factor when HQIV. |
| **fortran/camb.f90** | In `CAMB_ReadParams`: read `HQIV`, `hqiv_Ha_file`, `hqiv_cs2_fac`, `hqiv_l_cut`, `hqiv_beta`. |
| **inifiles/params_hqiv.ini** | New example inifile with `HQIV = T`, baryons-only, flat, and HQIV options. |

## New files (no patch to existing)

| File | Role |
|------|------|
| **camb_hqiv/patches/camb_hqiv.patch** | Single unified diff for all of the above (except new inifile). |
| **camb_hqiv/inifiles/params_hqiv.ini** | Example HQIV run (can be copied to CAMB `inifiles/`). |
| **camb_hqiv/run_hqiv_test.py** | Test script: age, first peak ℓ, low-ℓ suppression. |
| **camb_hqiv/hqiv_camb_wrapper.py** | Python wrapper: load modified camb, run HQIV, plot and save TT. |

## Optional (required if using Python CAMB from same repo)

| File | Role |
|------|------|
| **camb/model.py** | Add `HQIV`, `hqiv_Ha_file`, `hqiv_cs2_fac`, `hqiv_l_cut`, `hqiv_beta` to `CAMBparams._fields_` in the **same order** as in `fortran/model.f90` so the ctypes struct layout matches the Fortran type. See `patches/hqiv_python_model_snippet.py`. |

## Summary of physics changes (when HQIV = T)

1. **Background:** Pure baryons (`omch2=0`), flat; H(a) from `hqiv_Ha.txt`; dτ/da = 1/(a² H(a)); optional G(a) = G0 (Θ0/Θ(a))^0.6.
2. **Perturbations:** Same H(a) and optional G(a) in equations; optional effective sound-speed factor (placeholder ~0.95–1.05).
3. **Low-ℓ:** Exponential cutoff exp(-(ℓ/ℓ_cut)^1.8) for ℓ < ~80 (super-horizon damping).
4. **Recombination:** Unchanged (standard physics, no tired light).
