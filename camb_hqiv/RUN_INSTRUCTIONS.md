# HQIV CAMB — run instructions

## First principles (no DE language, age is prediction)

- **No separate dark-energy component.** Late-time acceleration comes entirely from the QI horizon motive force (Casimir-like cutoff on vacuum modes). The code reports:
  - *Effective late-time acceleration from QI horizon term: a₀ ~ 1.2×10⁻¹⁰ m/s² (matches observation)*
  - *Ω_total = 1 enforced by horizon equilibrium (no separate DE)*
- **Age is a prediction.** ~193 Gyr is the natural output of the QI axiom + varying G + horizon term. Do not tweak β or α_G to force 20 Gyr; β stays O(1) from the QI derivation (0.78). The older universe is a testable prediction (e.g. JWST galaxies get ~6–7× more proper time to form).

---

## 1. Background table (keep this)

Use the clean table with **H(a=1)/H₀ = 1** exactly:

```bash
cd sandbox && .venv/bin/python generate_hqiv_background.py
```

This writes `hqiv_Ha.txt`. Copy it to your CAMB run directory (or set `hqiv_Ha_file` in the .ini).

---

## 2. Quick test (<10 min on a typical PC)

In your params .ini (or Fortran `params_hqiv.ini`):

- `l_max_scalar = 500`
- `k_eta_max_scalar = 1200`
- `accuracy_boost = 0.8`

Then run CAMB (Fortran or Python). This should finish in under 10 minutes.

---

## 3. Full spectrum (l_max = 2500)

- Set `l_max_scalar = 2500`, `k_eta_max_scalar` as needed (e.g. 5000), `accuracy_boost = 1`.
- Conformal time η₀ is larger in HQIV (older universe) — the run takes longer. This is physics, not a bug.
- Run overnight or on a stronger machine.

---

## 4. Key numbers to output

When the run finishes, report:

- **Age today** (Gyr)
- **First-peak ℓ**
- **Low-ℓ suppression** (ℓ=2–30)
- **σ₈**
- **D(z=14) / D_ΛCDM**

Use `run_hqiv_covariant_test.py` (or your equivalent) to print these; the covariant extension uses the same clean table.
