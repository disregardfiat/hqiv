# HQIV Implementation in ECOSMOG

This directory contains the implementation plan and patches for testing HQIV (Horizon-Quantized Informational Vacuum) using the ECOSMOG N-body/PM simulation framework.

## Testable Predictions

### Core Predictions (from paper)

| Prediction | Description | Falsification Criterion |
|------------|-------------|------------------------|
| **Baryon fraction** | η = n_b/n_γ emerges from horizon statistics | Simulation must predict observed η |
| **CMB acoustic peaks** | Peak positions from modified growth | Peaks must align with Planck observations |
| **High-z proper time** | Extended proper time at high redshift | Must match JWST galaxy age constraints |
| **Temperature uniformity** | Horizon-scale mode coupling smooths CMB | σ_T/T consistent with CMB observations |
| **Structure growth** | Modified growth from horizon effects | Consistent with σ₈ and redshift-space distortions |

### Additional Predictions Testable with ECOSMOG

| Prediction | Description | Falsification Criterion | ECOSMOG Test |
|------------|-------------|------------------------|--------------|
| **Matter power spectrum** | Modified P(k) from horizon effects | P(k) shape matches observations | PM solver output |
| **Halo mass function** | Different high-z halo abundances | Match JWST galaxy cluster counts | Halo finder on N-body output |
| **Angular momentum excess** | Enhanced halo spin from vorticity | λ parameter higher than ΛCDM prediction | Compute J/E for halos |
| **Growth rate fσ₈** | Scale-dependent growth | Consistent with redshift-space distortions | Measure σ₈(z) and f(z) |
| **Void dynamics** | Different expansion in low-density regions | Void density profiles match observations | Identify voids, measure profiles |
| **Bullet Cluster test** | Lensing-gas separation from modified force | Recover observed separation | N-body + ray-tracing |
| **Correlation function** | Modified ξ(r) on horizon scales | Match SDSS/DESI correlation | Compute from particle positions |
| **Peculiar velocity field** | Different bulk flows | Consistent with velocity surveys | Measure v_bulk from simulation |

### Falsifiable Claims

1. **If CMB peaks are NOT consistent with Planck** → HQIV perturbation equations are wrong
2. **If high-z proper time is NOT extended** → Background expansion is wrong
3. **If halo mass function at high-z is NOT enhanced** → Structure formation prediction fails
4. **If Bullet Cluster separation is NOT reproduced** → Modified force law is incorrect
5. **If η does NOT emerge from early universe simulation** → Horizon statistics hypothesis fails
6. **If solar system tests show deviation** → High-acceleration limit is wrong (must match GR)

---

## ECOSMOG Background

ECOSMOG (Extended COmoving Symplectic MOnte-carlo Gadget) is a modified-gravity N-body code based on GADGET-2, with particle-mesh (PM) solvers for scale-dependent forces.

- **Repository**: https://github.com/MiXedCosmology/ECOSMOG
- **Key feature**: Already supports f(R), DGP, and general Horndeski gravity
- **Architecture**: PM solver on uniform grid + particle pusher

---

## Implementation Roadmap

### Phase 1: Background Module (1 week)

**File**: `src/background.c`

Replace standard Friedmann integration with HQIV H(a) table:

```c
// In background_init(), add:
if (HQIV_flag) {
    load_HQIV_table("hqiv_Ha.txt");
    // H(a) from table, not Friedmann equation
}

// Hubble function:
double H_of_a_HQIV(double a) {
    return interpolate_table(a, hqiv_a, hqiv_H_over_H0) * H0;
}
```

**Generate table (paper-aligned route)**:
```bash
cd ../../
python horizon_modes/python/bulk.py
# From the resulting lattice/background table, extract H(a) into ecosmog/hqiv_Ha.txt (a, H/H0)
```

### Phase 2: Modified Poisson Solver (2-4 weeks)

**File**: `src/poisson.c` or `src/forces.c`

The key modification: scale-dependent effective gravitational coupling.

```c
// Standard Poisson: ∇²Φ = 4πG ρ
// HQIV: ∇²Φ = 4πG_eff(a,k) ρ + horizon_term

double G_eff_HQIV(double a, double k) {
    double H_a = H_of_a_HQIV(a);
    double G_ratio = pow(H_a / H0, 0.6);  // α = 0.6 from QI
    
    // Horizon cutoff for k < 2π/Θ
    double Theta = 2 * c / H_a;  // in Mpc
    double k_horizon = 2 * M_PI / Theta;
    double horizon_factor = 1.0;
    
    if (k < k_horizon) {
        horizon_factor = exp(-k * Theta / (2 * M_PI));
    }
    
    return G0 * G_ratio * horizon_factor;
}
```

### Phase 3: Inertia Reduction in Particle Dynamics (2-4 weeks)

**File**: `src/predict.c` or `src/runge.c`

Modified inertia affects particle acceleration response:

```c
// In particle pusher, modify effective mass:
double inertia_factor(double a, double a_local) {
    // a_min from horizon scale
    double H_a = H_of_a_HQIV(a);
    double a_min = horizon_factor * c * H_a;
    
    if (fabs(a_local) < a_min) {
        return 0.35;  // floor to prevent instability
    }
    
    return 1.0 - a_min / fabs(a_local);
}

// Effective acceleration:
void accelerate_particle(int i, double a) {
    double a_local = compute_local_acceleration(i);  // from nearby particles
    double f_inertia = inertia_factor(a, a_local);
    
    // Standard: a = F/m
    // HQIV: a = F / (m * f_inertia)
    particle[i].acc = compute_force(i) / (particle[i].mass * f_inertia);
}
```

### Phase 4: Vorticity Source (4-8 weeks)

**File**: New file `src/vorticity.c`

The vorticity amplification term from the paper:

```c
// ∂ω/∂t = horizon_factor * H * (ω · ê_Θ)
// This is the most speculative part

void update_vorticity(double a, double dt) {
    double H_a = H_of_a_HQIV(a);
    
    for (int i = 0; i < N_particles; i++) {
        // Compute local vorticity from velocity curl
        vec3 omega = compute_vorticity(i);
        
        // Horizon direction (radial from observer)
        vec3 e_Theta = normalize(particle[i].position);
        
        // Amplification
        double amplification = horizon_factor * H_a * dot(omega, e_Theta);
        
        particle[i].vorticity = omega * (1.0 + amplification * dt);
    }
}
```

---

## Testing Strategy

### Test 1: Baryon Fraction from Early Universe

**Setup**: Start simulation at high temperature with only radiation and horizon-quantized modes.

**Procedure**:
1. Initialize radiation bath
2. Count modes that fit inside successive past light-cones
3. Track matter-antimatter asymmetry from horizon effects
4. Evolve to recombination

**Success criterion**: Output η consistent with observations

### Test 2: CMB Peak Structure

**Setup**: Run linear perturbation solver with HQIV modifications.

**Procedure**:
1. Generate H(a) table from background solver
2. Implement modified growth equations
3. Compute C_ℓ spectrum
4. Compare peak positions to Planck

**Success criterion**: Peaks consistent with Planck

### Test 3: High-z Age

**Setup**: Background integration only.

**Procedure**:
```python
# Implemented against the HQIV bulk/lattice background (see horizon_modes/python/bulk.py)
t_z = proper_time_at_z(a_arr, H_arr, z_target)
```

**Success criterion**: t(z) consistent with JWST galaxy age constraints

### Test 4: Temperature Uniformity

**Setup**: Track photon mode coupling across horizon scales.

**Procedure**:
1. Initialize CMB temperature field at recombination
2. Apply horizon-scale mode coupling
3. Measure temperature variance on different angular scales

**Success criterion**: σ_T/T consistent with CMB observations

---

## File Structure

```
ecosmog/
├── README.md                    # This file
├── hqiv_pm.py                   # Python PM simulator
├── params/
│   └── hqiv_params.ini          # ECOSMOG parameter file (future)
├── src/
│   ├── background_hqiv.c        # Modified background module (future)
│   ├── poisson_hqiv.c           # Modified Poisson solver (future)
│   ├── dynamics_hqiv.c          # Modified particle dynamics (future)
│   └── vorticity.c              # Vorticity source term (future)
├── tests/
│   └── (test scripts)
└── output/
    └── .gitkeep
```

## Current Status

**WORKING**: Python PM simulator (`hqiv_pm.py`) is functional and produces:
- Matter power spectrum P(k) evolution
- Structure growth tracking
- HQIV-specific modified gravity effects

### Key Predictions

1. **Extended universe age** - testable with stellar ages
2. **Modified structure growth** - testable with σ₈ measurements
3. **Higher H(z) at early times** - testable with BAO at high z
4. **Extended proper time at high-z** - allows older galaxies (JWST test)
5. **Enhanced halo spin** - testable with galaxy rotation curves
6. **Non-zero vorticity** in large-scale structure - testable with velocity field surveys

---

## Dependencies

- ECOSMOG base code: `git clone https://github.com/MiXedCosmology/ECOSMOG`
- FFTW3: `sudo apt install libfftw3-dev`
- GSL: `sudo apt install libgsl-dev`
- MPI: `sudo apt install libopenmpi-dev`
- Python: numpy, scipy, matplotlib

---

## Running the CLASS-HQIV code (Boltzmann / CMB)

CLASS with HQIV implements the HQIV background ($3H^2 - \gamma H = 8\pi G_{\rm eff}\,\rho$), inertia reduction, and vorticity source for CMB and linear structure. Use it for $C_\ell^{\rm TT}$, background tables, and age integration.

**The full CLASS and class_hqiv trees are not in this repository** (they are git-ignored). You must **download CLASS and apply our patches** before building.

### 1. Download CLASS and apply HQIV patches

Follow the instructions in **`class_hqiv_patches/README.md`** in this repository:

1. Clone CLASS: `git clone https://github.com/lesgourg/class_public.git`
2. Copy the HQIV source and header files from `class_hqiv_patches/source/` and `class_hqiv_patches/include/` into the CLASS `source/` and `include/` directories.
3. Add `hqiv.o` to the `SOURCE` line in the CLASS `Makefile`.
4. Copy the test `.ini` files from `class_hqiv_patches/` into the CLASS directory.

### 2. Build CLASS-HQIV

From the **CLASS** directory (after patching):

```bash
make clean
make class
```

Requirements: C and C++ compilers (gcc/g++), Make. Optional: HyRec, Recfast, Halofit (see CLASS docs and `external/`).

### 3. Run CLASS-HQIV

1. **Run with an HQIV parameter file** (e.g. `test_hqiv.ini` copied from `class_hqiv_patches/`):

   ```bash
   ./class test_hqiv.ini
   ```

   Output files use the prefix from `root` in the `.ini` (e.g. `root = hqiv_` → `hqiv_background.dat`).

2. **Key HQIV parameters** in the `.ini` file:

   - `HQIV = yes` — enable HQIV modifications
   - `HQIV_gamma = 0.4` — horizon term (≈ 0.35–0.45)
   - `HQIV_alpha = 0.60` — exponent in $G_{\rm eff}(a) = G_0(H/H_0)^\alpha$
   - `omega_cdm = 0.0` — baryon-only (no CDM)
   - `omega_lambda = 0.0` — no cosmological constant

3. **Optional: background only** (fast check): set `write_background = yes` in the `.ini` and inspect `*_background.dat`.

4. **Python (classy)**: `make classy` in the CLASS directory, then use `classy` from Python as usual.

---

## Build Instructions (ECOSMOG)

```bash
# Clone ECOSMOG
git clone https://github.com/MiXedCosmology/ECOSMOG.git
cd ECOSMOG

# Apply HQIV patches
cp ../ecosmog/src/*.c src/
cp ../ecosmog/params/hqiv_params.ini .

# Generate H(a) table (from HQIV bulk/lattice output; see horizon_modes/python/bulk.py)
cp ../hqiv_Ha.txt .

# Build
make clean && make

# Run test
mpirun -np 4 ./ECOSMOG hqiv_params.ini
```

---

## Expected Timeline

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| 1. Background | 1 week | H(a) table integration |
| 2. Poisson | 2-4 weeks | Modified force solver |
| 3. Dynamics | 2-4 weeks | Inertia reduction |
| 4. Vorticity | 4-8 weeks | Angular momentum evolution |
| 5. Testing | 2-4 weeks | All predictions validated |

**Total**: 3-6 months for complete implementation and testing

---

## Theoretical Background

The quantised inertia framework has seen renewed theoretical interest from the horizon physics community. Notably, Brodie (2026) derives a MOND-like modification to inertia from Jacobson thermodynamics, finding the critical acceleration scale within 9% of Milgrom's empirical value. This independent theoretical derivation from thermodynamic first principles connects the QI/horizon physics programme to the broader MOND literature and strengthens the case for horizon-based modifications to inertia.

## References

1. ECOSMOG paper: Li et al., JCAP 2013, "ECOSMOG: An Efficient Code for Simulating Modified Gravity"
2. GADGET-2 manual: Springel, MNRAS 2005
3. QI framework: McCulloch, arXiv:1610.06787 and subsequent works
4. Thermodynamic derivation: Keith Brodie, Zenodo:10.5281/zenodo.18706746 (2026) — derives MOND-like acceleration from Jacobson thermodynamics
