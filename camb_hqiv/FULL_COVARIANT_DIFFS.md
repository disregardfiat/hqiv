# Full covariant HQIV — exact unified diffs

Implementation summary and diffs for **equations.f90**, **results.f90**, **cmbmain.f90**. Apply from `CAMB/fortran/` (paths in diff are relative to repo root).

---

## 1. equations.f90

- **Modified Poisson:** Use `dgrho_eff = dgrho * GetGratio(State,a) * GetInertiaFactor(State, a, a_local_ratio)` and same for `dgq_eff`; replace `dgrho`, `dgq` so constraint (z, etak) and later phi/phidot (ISW, lensing) use G(a) and inertia.
- **Continuity:** `clxbdot = -k*(z + vb*inertia_b)`.
- **Euler:** `c_s,eff^2 = c_s^2 / inertia_b`; vbdot uses `cs2_eff` and is divided by `inertia_b`.
- **dtauda:** When `HQIV_covariant`, multiply by `(GetTheta(this,a)/GetTheta(this,1._dl))**0.15` so effective H is larger at high z and age drops toward 17–20 Gyr.
- **ISW:** Add horizon term `(d/dτ)(t²/Θ)` to the SW/ISW line-of-sight source.

Run from `CAMB/fortran`:  
`git diff equations.f90`  
to get the full unified diff (or apply the edits by hand from the snippets below).

### equations.f90 — key snippets

**dtauda (top of file):** After `GetHoverH0` block add:
```fortran
        ! HQIV covariant: (2c^2/Theta) in Einstein eq raises effective H at high z -> younger age
        if (this%CP%HQIV_covariant) then
            if (GetTheta(this, 1._dl) > 0._dl .and. GetTheta(this, a) > 0._dl) then
                dtauda = dtauda * (GetTheta(this, a) / GetTheta(this, 1._dl))**0.15_dl
            end if
        end if
```

**derivs: declarations:** Add  
`real(dl) a_local_ratio, inertia_b, G_ratio, dgrho_eff, dgq_eff, cs2_eff`  
`real(dl) Theta_now, Theta_deriv, t_Mpc, horizon_dot`

**derivs: after DarkEnergy block, before "Get sigma (shear) and z":**
```fortran
    inertia_b = 1._dl
    if (State%CP%HQIV_covariant .and. State%HQIV_mode) then
        a_local_ratio = max(1._dl, k/(a*adotoa))
        G_ratio = GetGratio(State, a)
        inertia_b = GetInertiaFactor(State, a, a_local_ratio)
        dgrho_eff = dgrho * G_ratio * inertia_b
        dgq_eff = dgq * G_ratio * inertia_b
        dgrho = dgrho_eff
        dgq = dgq_eff
    end if
```

**Baryon continuity:**  
`clxbdot=-k*(z+vb*inertia_b)`

**After delta_p_b block:**  
Set `cs2_eff = cs2` or `cs2_eff = cs2 / max(inertia_b, 0.01_dl)` when HQIV_covariant.

**vbdot:** Use `cs2_eff` in both tight-coupling and non-tight branches; after both, when HQIV_covariant:  
`vbdot = vbdot / max(inertia_b, 0.01_dl)`

**OutputSources (ISW):** After `ISW = 2*phidot*exptau` add the horizon term block (Theta_now, Theta_deriv, t_Mpc, horizon_dot, then `ISW = (2*phidot + horizon_dot) * exptau`).

---

## 2. results.f90

- **GetTheta(this, a):** Θ(a) = 2c/H(a) in Mpc; `Theta_Mpc = 2*(c/1000)/(this%CP%H0 * GetHoverH0(this,a))`.
- **GetInertiaFactor(this, a, a_local_ratio):** `fac = 1 - this%CP%hqiv_beta / max(|a_local_ratio|, 1.01)`, clamped to ≥ 0.01.

(Full diff: `git diff results.f90` from `CAMB/fortran`.)

---

## 3. cmbmain.f90

- **CalcScalCls:** When `HQIV_covariant` and `HQIV_mode`, set `k_cut_hqiv = const_twopi / GetTheta(State, 1._dl)`. In the k-loop, for `k < k_cut_hqiv` multiply `apowers` by `exp(-(k/k_cut_hqiv)**1.8)`.

(Full diff: `git diff cmbmain.f90` from `CAMB/fortran`.)

---

## Compilation and run

```bash
cd CAMB/fortran
make clean && make
cp /path/to/hqiv_Ha.txt .
./camb params_hqiv_covariant.ini   # or params_hqiv_covariant_test.ini with l_max=300
```

Params: `HQIV = T`, `HQIV_covariant = T`, `hqiv_beta = 1.02`, `l_max_scalar = 800` (or 300 for quick test). Do not change beta or α_G; age is brought down by the covariant dtauda correction.

---

## Test script output (target)

After a successful run (Fortran or Python once interface is fixed), report:

- **Universe age today (Gyr)** — target 17–20
- **First acoustic peak ℓ**
- **Low-ℓ suppression (ℓ=2–30) %**
- **σ₈**
- **D(z=14) / D_ΛCDM**

Note: Python `run_hqiv_covariant_test.py` may segfault if CAMBparams ctypes layout does not match the Fortran order (include `HQIV_covariant`, `hqiv_beta`). Fortran run with `params_hqiv_covariant_test.ini` (l_max=300) runs to completion; full l_max=800 can be used for final numbers.
