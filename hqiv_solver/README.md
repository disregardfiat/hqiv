# HQIV Perturbation Solver

A pure-Python (NumPy + SciPy) Fourier-space integrator for the Horizon-Quantized Informational Vacuum (HQIV) framework. This solver evolves the exact scalar and vector perturbation equations from deep radiation era (a ≈ 10⁻⁹) through recombination to today.

## Paper Reference

**"Horizon-Quantized Informational Vacuum (HQIV): A Covariant Baryon-Only Cosmological Framework from Quantised Inertia"**

This solver implements the key equations from the paper:
- **Background dynamics** (Section 5, Eqs. 13-15)
- **Scalar perturbations** (Section 7, Eqs. 9-11)
- **Vector perturbations** (Section 7, Eq. 12)
- **Observables** (Section 8)

## Features

- ✅ Background cosmology with HQIV horizon term
- ✅ Scalar perturbations with modified inertia
- ✅ Vector (vorticity) perturbations with horizon coupling
- ✅ CMB angular power spectrum C_ℓ^TT up to ℓ = 500
- ✅ Growth factor and σ₈ computation
- ✅ Toggle between HQIV and standard ΛCDM

## Installation

```bash
# Clone the repository
git clone https://github.com/disregardfiat/hqiv.git
cd hqiv

# Install dependencies
pip install numpy scipy matplotlib

# The solver is in hqiv_solver/
```

## Quick Start

```python
from hqiv_solver import HQIVPerturbations

# Create HQIV solver with best-fit parameters
solver = HQIVPerturbations(
    H0=73.2,           # Hubble constant [km/s/Mpc]
    Om_m=0.031,        # Matter density (baryons only)
    hqiv_on=True,      # Enable HQIV modifications
    beta=1.02,         # Horizon-smoothing parameter
    Om_h=1.00,         # Horizon density parameter
    n_h=1.04           # Horizon dilution exponent
)

# Run the solver
solver.run()

# Compare with ΛCDM
solver.compare_with_lcdm()

# Plot results
solver.plot_results()

# Print summary
solver.summary()
```

## Key Predictions

| Observable | HQIV Prediction | ΛCDM | Notes |
|------------|-----------------|------|-------|
| Universe age | ~17 Gyr | 13.8 Gyr | Addresses JWST timing |
| H₀ | 73.2 km/s/Mpc | 67.4 km/s/Mpc | Matches SH0ES |
| Growth factor | ~0.36× ΛCDM | 1.0 | Falsifiable prediction |
| Vorticity | Amplified (β > -2) | Decays as a⁻² | Unique signature |
| CMB peak | Shifted | ℓ ~ 220-280 | Tests acoustic physics |

## Module Structure

```
hqiv_solver/
├── __init__.py          # Package initialization
├── background.py        # Background cosmology (H(a), η(a), G_eff(a))
├── scalars.py           # Scalar perturbations (δ, θ, Φ)
├── vectors.py           # Vector perturbations (vorticity ω)
├── observables.py       # CMB spectrum, P(k), σ₈
├── hqiv_solver.py       # Main HQIVPerturbations class
├── README.md            # This file
└── notebooks/
    └── demo.ipynb       # Example Jupyter notebook
```

## Equations Implemented

### Background (Paper Eq. 15)
```
H² = H₀² [Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_h a⁻ⁿ]
```

### Scalar Perturbations (Paper Eqs. 9-11)

**Continuity:**
```
δ' + θ = -3 Φ'
```

**Euler (modified inertia):**
```
θ' + ℋ θ = -k² Φ / f(α, φ)
```
where `f(α, φ) = √(1 - cφ/α)` (HQIV form)

**Poisson (with horizon correction):**
```
k² Φ = -(3/2) ℋ² Ω_m(a) δ a² + horizon_term × a²
```

### Vector Perturbations (Paper Eq. 12)

**Vorticity evolution:**
```
∂ω/∂t + (v·∇)ω = β H (ω · ê_Θ)
```

In conformal time:
```
ω' + ℋ ω = β(a) a H(a) (ω · ê_Θ)
```

## Validation

### ΛCDM Limit
When `hqiv_on=False`, the solver recovers standard ΛCDM:
- Growth factor D(a) ∝ a in matter domination
- Vorticity decays as ω ∝ a⁻²
- CMB peak at ℓ ~ 220-280

### HQIV Mode
When `hqiv_on=True`:
- Growth suppressed by factor ~0.36
- Vorticity can grow or remain constant
- CMB peak position may shift

## Example Output

```
============================================================
HQIV Perturbation Solver - Summary
============================================================

Parameters:
  H0: 73.2
  Om_m: 0.031
  hqiv_on: True
  beta: 1.02
  Om_h: 1.00
  n_h: 1.04

Background:
  Universe age: 17.1 Gyr
  Sound horizon: 145.2 Mpc

Scalar Perturbations:
  Growth factor at a=0.5: 0.4521

Vector Perturbations:
  Vorticity growth exponent: -0.85
  Final amplification: 3.45e+02

Observables:
  CMB first peak: ℓ = 235

Comparison with ΛCDM:
  Growth suppression at a=0.5: 0.38
  CMB peak shift: Δℓ = 15
============================================================
```

## Dependencies

- Python 3.11+
- NumPy
- SciPy
- Matplotlib (for plotting)

## License

MIT License - See LICENSE file for details.

## Citation

If you use this code, please cite the accompanying paper:

```bibtex
@article{ettinger2026hqiv,
  title={Horizon-Quantized Informational Vacuum (HQIV): 
         A Covariant Baryon-Only Cosmological Framework 
         from Quantised Inertia},
  author={Ettinger Jr, Steven},
  journal={arXiv preprint},
  year={2026}
}
```

## References

1. McCulloch, M.E. (2007). "Modelling the Pioneer anomaly as modified inertia"
2. McCulloch, M.E. (2016). "Quantised inertia from relativity and the uncertainty principle"
3. Brodie, K. (2026). "Thermodynamic derivation of MOND from horizon physics"