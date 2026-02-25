# σ₈ and full CLASS output (peak-scan best-fit input)

## What is σ₈ (sigma8)?

**σ₈** is the **RMS amplitude of matter density fluctuations** in spheres of radius **8 Mpc/h** at **z = 0**, in the linear regime. It is defined as

- σ₈² = ∫₀^∞ (3 j₁(kR)/(kR))² P(k) k²/(2π²) dk  with R = 8 Mpc/h,

where P(k) is the linear matter power spectrum and j₁ is the first spherical Bessel function. It is a standard observable used to normalize the amplitude of structure formation (e.g. Planck reports σ₈ ≈ 0.81 in ΛCDM; weak lensing and cluster counts often quote σ₈ or S₈ = σ₈ √(Ω_m/0.3)).

---

## Input used (peak-alignment scan best-fit)

| Parameter    | Value  |
|-------------|--------|
| gamma (HQIV) | 0.36  |
| omega_b     | 0.025 |
| h           | 0.732 |
| HQIV_alpha  | 0.60  |
| omega_cdm   | 0     |
| A_s         | 2.1e-9 |
| n_s         | 0.96  |
| tau_reio    | 0.054 |

---

## Key outputs

| Quantity | Value |
|----------|--------|
| **σ₈** (total matter) | **0.7629** (computed till k = 1.57 h/Mpc) |
| Age (Gyr) | 19.987 |
| Conformal age (Mpc) | 26315.7 |
| Ω_Λ (effective closure) | 0.833 |
| Ω_b | 0.0467 |
| z_eq | 596.7 |
| Sound horizon at rec (Mpc) | 192.75 |
| 100 θ_s | 0.743 |
| 100 θ_* | 0.747 |

Full run log is in `full_output2.txt`.
