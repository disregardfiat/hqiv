# Full covariant HQIV extension — files to modify

All paths relative to the **CAMB repository root**. This extends the existing HQIV fork (tabulated H(a), Ω_cdm=0, G(a), low-ℓ damping) with the **covariant HQIV** formulation.

## Philosophy

- **HQIV_covariant = F**: code behaves exactly like standard CAMB (and like existing HQIV when only `HQIV = T`).
- **HQIV_covariant = T**: apply horizon-modified metric, modified Einstein term (2c²/Θ)g_μν, inertia reduction m_i = m_g[1 − a_min/|a_local|], effective sound speed, varying G in Poisson/growth, and horizon cutoff in SW/ISW source terms.
- Every change commented: `! HQIV covariant extension — Steven Ettinger + Mr 4.20`

---

## 1. fortran/model.f90

- Add to `CAMBparams`:
  - `logical :: HQIV_covariant = .false.`
- Ensure existing: `HQIV`, `hqiv_Ha_file`, `hqiv_cs2_fac`, `hqiv_l_cut`, `hqiv_beta` (β ≈ 0.78 for a_min = β c H).

---

## 2. fortran/results.f90

- Add:
  - `GetTheta(this, a)`: Θ(a) = 2c/H(a) from tabulated H (for horizon term and cutoff scale).
  - `GetInertiaFactor(this, a, a_local)` or equivalent: inertia reduction factor for continuity/Euler.
  - When `HQIV_covariant = T`, ensure σ8 and growth use G(a) where appropriate (varying G in growth/transfer).
- Optional: store Θ(a) on the same spline as H(a) for fast evaluation in source terms.

---

## 3. fortran/equations.f90

- **Background / dtauda**: already uses tabulated H(a) when HQIV; no change for covariant (same background table).
- **Modified Einstein term**: the term (2c²/Θ)g_μν enters as an effective contribution in the 00 and ij equations. In the perturbation equations this implies:
  - An effective pressure/density contribution from the horizon (constant in background; can be folded into the existing effective DE from the table).
  - In the Poisson equation (relating Φ to δρ): include varying G(a) and, if needed, the linearised horizon term.
- **Inertia reduction**: in the baryon-photon fluid:
  - Effective sound speed: `c_s,eff² = c_s² / (1 − β c H / |a_fluid|)` with |a_fluid| a safe proxy (e.g. max(c_s k/a, c H)).
  - Modified continuity and Euler: multiply inertial terms by factor `[1 − a_min/|a_local|]` (i.e. m_i/m_g) where a_min = β c H(t).
- **Horizon cutoff**: in the CMB source terms (Sachs-Wolfe + ISW), apply cutoff for modes with k < 2π/Θ(τ) (or k < k_cut(τ) = 2π/Θ(τ)) in the integrand.

---

## 4. fortran/camb.f90

- In `CAMB_ReadParams`: read `HQIV_covariant` (default F). When T, require `HQIV = T` and same `hqiv_Ha_file`.

---

## 5. fortran/subroutines.f90 (or equivalent: source integration)

- Where the C_l source functions are integrated over τ/k: apply the horizon cutoff factor (e.g. exp(-(k/k_cut(τ))^1.8) or step) for k < 2π/Θ(τ) when `HQIV_covariant = T`.
- Alternatively, the cutoff can be applied in the same place as the existing low-ℓ damping (in results or in the Cl assembly) using k_cut = 2π/Θ(τ_dec) for the relevant τ.

---

## 6. inifiles/params_hqiv.ini (or new params_hqiv_covariant.ini)

- Add line: `HQIV_covariant = T` when running the full covariant case.
- Keep `HQIV = T`, `hqiv_Ha_file = hqiv_Ha.txt`, and other HQIV options.

---

## 7. Python: camb/model.py (if using Python interface)

- Add `HQIV_covariant` (c_bool) to `CAMBparams._fields_` in the **same order** as in `fortran/model.f90`.

---

## 8. New files (no modification to existing)

| File | Role |
|------|------|
| **patches_covariant/01_model.f90.patch** | Unified diff: add `HQIV_covariant` to CAMBparams. |
| **patches_covariant/02_camb.f90.patch** | Unified diff: read `HQIV_covariant` in CAMB_ReadParams. |
| **patches_covariant/03_results.f90.patch** | Unified diff: comment line; paste GetTheta/GetInertiaFactor from snippet. |
| **patches_covariant/04_equations.f90.patch** | Unified diff: comment line; paste metric/inertia/c_s_eff/cutoff from snippet. |
| **patches_covariant/results_covariant_snippet.f90** | GetTheta(a), GetInertiaFactor(a, a_local_scale) to paste in results.f90. |
| **patches_covariant/equations_covariant_snippet.f90** | Where to add c_s,eff, inertia, G(a), horizon cutoff in equations.f90. |
| **patches_covariant/source_horizon_cutoff_snippet.f90** | Horizon cutoff k < 2π/Θ in CMB source integrand. |
| **patches_covariant/python_model_covariant_snippet.py** | Add `HQIV_covariant` to Python CAMBparams._fields_. |
| **run_hqiv_covariant_test.py** | Python test: age, first peak ℓ, low-ℓ suppression, σ8, D(z=14) vs ΛCDM. |
| **inifiles/params_hqiv_covariant.ini** | Example .ini with HQIV_covariant = T and get_transfer for σ8/D(z). |
| **COMPILE_AND_RUN.md** | Compilation and run instructions. |

---

## Summary of covariant physics (when HQIV_covariant = T)

1. **Metric**: Horizon-modified term (c² t²/Θ_local) in g_00 to 1PN; Θ_local = 2c/|a_local|.
2. **Einstein**: G_μν + (2c²/Θ)g_μν = (8πG/c⁴)T_μν; Θ = 2c/H in background.
3. **Inertia**: m_i = m_g[1 − a_min/|a_local|], a_min = β c H; applied in continuity/Euler and in c_s,eff².
4. **Varying G**: G(a) = G0 (Θ0/Θ(a))^0.6 in Poisson and growth.
5. **Horizon cutoff**: k < 2π/Θ in SW/ISW source terms (or equivalent low-k damping).
6. **Outputs**: TT, TE, EE, lensed, matter P(k), σ8, growth D(z) — all comparable to Planck.
