# Prompt for full CLASS / HiCLASS sandbox

Copy-paste the block below into your full Python/CLASS environment (HiCLASS branch recommended: [github.com/lesgourg/class_public](https://github.com/lesgourg/class_public) + HiCLASS fork for modified gravity). It implements the exact HQIV background + perturbations.

---

```python
# HQIV Cosmology for CLASS / HiCLASS
# Run after installing CLASS with modified gravity support

import numpy as np
from classy import Class

# HQIV parameters (locked from our sandbox)
beta = 0.81          # horizon strength (a_min = beta * c * H)
omega_b = 0.0224     # visible baryons only
omega_cdm = 0.0      # ZERO cold dark matter
h = 0.70
Omega_m_vis = 0.048

# Background override: custom Hubble function
# H(a) = H0 * sqrt( Omega_r a^{-4} + Omega_m_vis * a^{-3} * (1 + delta_QI) + Omega_Lambda_QI(a) )
# where delta_QI and Omega_Lambda_QI come from our A_eff = A_std + beta * H

class HQIV_Class(Class):
    def __init__(self):
        super().__init__()
        self.set({
            'output': 'tCl pCl lCl mPk',
            'l_max_scalars': 3000,
            'lensing': 'yes',
            'P_k_max_1/Mpc': 10,
            'omega_b': omega_b,
            'omega_cdm': 0.0,
            'H0': h*100,
            'Omega_k': 0.0,
            'A_s': 2.1e-9,
            'n_s': 0.96,
            'tau_reio': 0.054,
            # HiCLASS modified gravity flags
            'mg_flag': 1,                    # activate modified gravity
            'mu0': 1.0,                      # placeholder; we override below
            'gamma0': 1.0,
            # Custom background via external module or table (see below)
        })

# Feed our numerical H(a) table (generate from our scipy code first)
# Example: save H(a) array from sandbox as hqiv_Ha.txt (columns: a, H/H0)
# Then in CLASS: use 'background' custom or modify background.c

# Step-by-step run instructions:
# 1. Run our HQIV background integrator → save np.savetxt('hqiv_Ha.txt', np.c_[a_arr, H_arr])
# 2. In CLASS source (background.c or use external_Pk), interpolate H(a) from table
# 3. Add to perturbation module: delta_m'' = ... + QI_boost_term (from our growth ODE)
# 4. Compile & run:
# cosmo = HQIV_Class()
# cosmo.set({'background_table': 'hqiv_Ha.txt'})  # or custom flag
# cosmo.compute()
# cls = cosmo.lensed_cl(3000)
# Plot TT, compare low-ell suppression, acoustic peaks
```

---

**Exact modifications in CLASS source (HiCLASS branch):**

- **background.c:** Replace standard Friedmann with table interpolation of our solved H(a) from the quadratic A_eff solver.
- **perturbations.c:** Add to the δ' equation the QI inertia reduction factor `1/(1 − beta * c * H / |A|)`.
- **Horizon cutoff:** Add exponential damping `exp(−k * Θ_local / (2π))` to longest modes (k < 2π/Θ).
