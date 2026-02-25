# Full run at cost minimum (Omega_m = 0.0191)

Deep run at the multipole-delta cost minimum from `background_cost_scan.py` (0.015–0.03 rescan). Full outputs: Cℓ, P(k), background, thermodynamics, perturbations.

## Input

| Parameter | Value |
|-----------|--------|
| Omega_m   | 0.0191 (baryon-only: omega_b = Omega_m × h²) |
| omega_b   | 0.01017259 |
| h         | 0.73 |
| T_cmb     | 2.725 K |
| HQIV      | yes (gamma=0.4, alpha=0.6, chi=0.172, fmin=0.01) |
| gauge     | newtonian |

## Output files

| File | Description |
|------|-------------|
| `run.ini` | Input parameters |
| `full_run.log` | Full CLASS stdout (budget, diagnostics, sigma8) |
| `hqiv_min*_cl.dat` | CMB Cℓ (TT; l_max=2500) |
| `hqiv_min*_pk.dat` | Matter P(k) at z_pk=0 |
| `hqiv_min*_background.dat` | Background vs z (H, Omega's, etc.) |
| `hqiv_min*_thermodynamics.dat` | Visibility, sound horizon, theta_s, etc. |

(Number suffix 00/01 from CLASS output numbering.)

## Full picture at 0.0191 (from full_run.log)

### Time and expansion
- **Age (wall-clock)** = 51.18 Gyr
- **Conformal age** = 49597 Mpc
- **H_input (h=0.73)** = 73.0 km/s/Mpc
- **H_actual(z=0)** = 16.1 km/s/Mpc → ratio H_actual/H_input = 0.220
- **1/H_actual(z=0)** = 60.79 Gyr (Hubble time)
- **Time dilation:** wall-clock age / ΛCDM age(H_input) ≈ 3.96× (universe ~4× older than ΛCDM with same h says)

### Energy budget at z=0
- **Omega_m** = 0.393 (physical; large because H_actual ≪ H_input)
- **Omega_r** = 1.62e-3
- **Omega_k** = 0
- **Omega_horizon (γH)** = 0.605 (replaces dark energy)
- **Sum** = 1.0

### HQIV diagnostics
- **rho_horizon / ΔE_rad (z_rec→0)** = 0.34 (horizon captured fraction)
- **E(z=0)/E(z_rec)** (comoving) = 0.46
- **z_eq** (radiation/matter) = 242.4

### CMB / thermo
- **z_rec** = 1093.6, **rs_rec** = 218.13 Mpc (comoving sound horizon)
- **D_A(z_rec)** = 44.94 Mpc
- **100×theta_s** = 0.4434, **100×theta_*** = 0.4452
- **z_drag** = 1004.6, **rs_drag** = 235.28 Mpc
- **z_reio** = 5.33

### Linear spectra
- **sigma8** = 0.0987 (total matter, to k = 1.65 h/Mpc)

## Re-run

```bash
/path/to/class_public/class run.ini
```
Or: `CLASS_EXEC=/path/to/class python3 ../background_cost_scan.py ...` for scans.
