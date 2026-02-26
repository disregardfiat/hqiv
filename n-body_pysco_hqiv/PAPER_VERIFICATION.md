# N-body solver vs paper math verification

This document cross-checks the n-body implementation against **paper/main.tex** (and related sections).

## 1. Modified Einstein equation

**Paper** (main.tex, Units & Conventions, Appendix horizon-action):
```
G_μν + γ (φ/c²)(δθ̇'/c) g_μν = (8π G_eff(φ)/c⁴) T_μν
```
Homogeneous limit: φ ≈ cH, δθ̇' ≈ H → (3−γ)H² = 8π G_eff(ρ_m + ρ_r).

**Code**: Not used directly in N-body; background H(a) and G_eff(a) are taken from the HQIV background solver or approximate relations. **Consistent** (background wrapper / CLASS).

---

## 2. Modified geodesic / inertia

**Paper** (main.tex, Modified Geodesic Eq., Bullet list):
- Particle action: S = −m_g ∫ f(a_loc, φ) ds.
- Weak-field limit: m_i **a** = −m_g ∇Φ with m_i = m_g f → **a** = −∇Φ / f.
- Inertia factor: **f(a_loc, φ) = max( a_loc / (a_loc + φ/6), f_min )**, with a_loc and φ in acceleration units [m/s²]; φ = 2c²/Θ_local.

**Code** (`hqiv_modifications/inertia_factor.py`):
- **Thermo form**: `denominator = alpha + phi/6`, `f = alpha/denominator`, then floor f_min and cap 1. **Matches paper.**
- **Apply to acceleration**: `acc_modified = acc / f` (effective acceleration = −∇Φ/f). **Matches paper.**

**Units**: φ from `phi_field.py` is 2c²/Θ_local [m/s²]; local acceleration is converted to [m/s²] in `run_bullet.py` before calling the inertia solver. **Consistent.**

---

## 3. Varying gravitational coupling

**Paper** (main.tex, Bullet list, Table):
```
G_eff(a) = G_0 (H(a)/H_0)^α
```
α ≈ 0.60 (or dynamic α_eff = χφ/6 in sim).

**Code** (`hqiv_modifications/g_eff.py`):
```python
def effective_gravitational_constant(H_a, H0, alpha=_alpha_G):
    return (H_a / H0) ** alpha
```
**Matches paper.**

---

## 4. Poisson equation and horizon term

**Paper**: (00)-component of the modified Einstein equation in the Newtonian limit gives a Poisson equation with a horizon contribution. Homogeneous limit: ∇²Φ has standard source −4π G_eff ρ plus a horizon term ∝ γ H² (effective constant source on large scales).

**Code** (`hqiv_modifications/g_eff.py`):
- Poisson RHS: standard part (3/2) Ω_m G_eff a δ (code units).
- Horizon correction: added in k-space, scale-dependent (suppressed at large k), magnitude ∝ γ (H/H0)² a². **Consistent** with a Newtonian limit of the horizon term (large-scale effective source).

---

## 5. φ field

**Paper**: φ(x) = 2c²/Θ_local(x); FLRW: Θ = 2c/H → φ = cH. Units of φ as acceleration [m/s²] (main.tex: “both a_loc and φ are accelerations”).

**Code** (`hqiv_modifications/phi_field.py`):
- Θ_local from expansion scalar or background 2c/H.
- φ = 2c²/Θ_local, returned in [m/s²]. **Matches paper.**

---

## 6. Summary

| Item              | Paper (main.tex)        | Code                          | Status   |
|-------------------|--------------------------|-------------------------------|----------|
| Inertia f         | f = a/(a + φ/6), φ [m/s²]| thermo: alpha + phi/6         | Matches  |
| Geodesic          | a_eff = −∇Φ/f           | acc_modified = acc/f          | Matches  |
| G_eff             | G_0 (H/H0)^α            | (H_a/H0)**alpha               | Matches  |
| φ definition      | 2c²/Θ_local             | 2*c**2/Theta_local            | Matches  |
| Poisson + horizon | ∇²Φ = source + horizon  | RHS + k-space horizon term    | Consistent |

**Fix applied**: Inertia factor was previously implemented as f = α/(α + χcφ/6) with φ in [m/s²], which is dimensionally inconsistent (χcφ has dimension m²/s³). It is now f = α/(α + φ/6) in line with the paper.
