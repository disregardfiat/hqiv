# Full covariant HQIV — files to modify

Implementation of Horizon-Quantized Informational Vacuum (McCulloch Quantised Inertia) with **dynamic horizon term only** and full covariant extensions.

---

## CRITICAL: Dynamic horizon term only

- **A_eff = A_std + beta × H(t)²** (not H₀²)
- **beta = 1.02** (fixed O(1) from QI axiom)
- Gives universe age **~17–20 Gyr**

---

## 1. Background (Python)

| File | Role |
|------|------|
| **horizon_modes/python/bulk.py** | Main HQIV background + lattice generator used in the paper; from its output you can extract H(a) and write `hqiv_Ha.txt` (a, H_over_H0) for CAMB. (Legacy: `sandbox/generate_hqiv_background.py` produced an earlier toy-background `hqiv_Ha.txt`.) |

---

## 2. CAMB Fortran — params and I/O

| File | Modifications |
|------|---------------|
| **fortran/model.f90** | Default `hqiv_beta = 1.02_dl`. Keep `HQIV`, `hqiv_Ha_file`, `hqiv_cs2_fac`, `hqiv_l_cut`, `HQIV_covariant`. |
| **fortran/camb.f90** | Read `HQIV_covariant`, `hqiv_beta` (already present). Ensure `hqiv_beta` read when HQIV=T. |

---

## 3. CAMB Fortran — background and Θ(a), inertia

| File | Modifications |
|------|---------------|
| **fortran/results.f90** | `GetTheta(this, a)`: Θ(a)=2c/H(a) from tabulated H. `GetInertiaFactor(this, a, a_local)`: m_i/m_g = 1 − beta×c×H/(\|a_local\|). When HQIV_covariant=T, use in growth/σ8. |
| **fortran/equations.f90** | dtauda: unchanged (uses table). When HQIV_covariant: modified Poisson (G(a)), inertia in continuity/Euler, c_s,eff² = c_s²/(1 − beta×c×H/\|a_fluid\|), horizon term in 00/ij. |

---

## 4. CAMB Fortran — horizon cutoff in CMB sources

| File | Modifications |
|------|---------------|
| **fortran/cmbmain.f90** (or **equations.f90** source section) | Apply super-horizon cutoff: **exp(-(k/k_cut)^1.8)** for k < k_cut = 2π/Θ₀ (low-ℓ modes). Apply when building source or when forming Δ_p_l_k. |

---

## 5. Inifiles and Python

| File | Role |
|------|------|
| **camb_hqiv/inifiles/params_hqiv_covariant.ini** | `HQIV = T`, `HQIV_covariant = T`, `hqiv_beta = 1.02`, `hqiv_Ha_file = hqiv_Ha.txt`. Get transfer for σ8, D(z). |
| **CAMB/camb/model.py** (if Python interface) | Add `HQIV_covariant` (and ensure `hqiv_beta`) in CAMBparams `_fields_` in same order as Fortran. |
| **camb_hqiv/run_hqiv_covariant_test.py** | Run with hqiv_Ha.txt; print age (~17–20 Gyr), first-peak ℓ, low-ℓ suppression, σ8, D(z=14)/D_ΛCDM. |

---

## 6. Summary of physics (HQIV_covariant = T)

1. **Background**: Tabulated H(a) from hqiv_Ha.txt (generated with **A_eff = A_std + beta×H(t)²**, beta=1.02). H(a=1)=H₀. Varying G(a)=G₀(Θ₀/Θ(a))^0.6.
2. **Metric**: ds² = −(1+2Φ+c²t²/Θ_local)c²dt² + a²(1−2Φ)δ_ij dx^i dx^j; Θ_local = 2c/|a_local|.
3. **Einstein**: G_μν + (2c²/Θ)g_μν = (8πG/c⁴)T_μν.
4. **Inertia**: m_i = m_g[1 − beta×c×H/|a_local|]; in fluid: c_s,eff² = c_s²/(1 − beta×c×H/|a_fluid|).
5. **Perturbations**: Modified continuity, Euler, Poisson with G(a) and inertia factor.
6. **Horizon cutoff**: exp(-(k/k_cut)^1.8), k_cut = 2π/Θ₀, for low-ℓ.

---

## Order of patches

1. **Background generator** (Python) — dynamic term, beta=1.02  
2. **model.f90** — hqiv_beta default 1.02  
3. **camb.f90** — confirm hqiv_beta read  
4. **results.f90** — GetTheta, GetInertiaFactor  
5. **equations.f90** — covariant terms (G(a), inertia, c_s,eff, horizon)  
6. **cmbmain.f90** or **equations.f90** — horizon cutoff in source/transfer  
7. **params_hqiv_covariant.ini** — HQIV_covariant=T, hqiv_beta=1.02  
8. **run_hqiv_covariant_test.py** — age, peak, low-ℓ, σ8, D(z=14)  
9. **COMPILE_AND_RUN.md** — build and test
