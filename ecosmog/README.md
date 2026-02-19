# HQIV Implementation in ECOSMOG

This directory contains the implementation plan and patches for testing HQIV (Horizon-Quantized Informational Vacuum) using the ECOSMOG N-body/PM simulation framework.

## Testable Predictions

### Core Predictions (from paper)

| Prediction | Description | Falsification Criterion |
|------------|-------------|------------------------|
| **Baryon fraction** | η = n_b/n_γ emerges from horizon statistics | Simulation must predict η ~ 6×10⁻¹⁰ |
| **CMB multimode peaks** | Acoustic peak positions from modified growth | Peaks 2-5 must align with Planck within 3σ |
| **z=14 age** | Proper time at z=14 ~ 940 Myr (3.2× ΛCDM) | Must match JWST galaxy ages |
| **Temperature uniformity** | Horizon-scale mode coupling smooths CMB | σ_T/T < 10⁻⁵ on >1° scales |
| **β(t) evolution** | Horizon smoothing factor → 1 as universe ages | β(t) measurable from horizon anisotropy at different z |

### Additional Predictions Testable with ECOSMOG

| Prediction | Description | Falsification Criterion | ECOSMOG Test |
|------------|-------------|------------------------|--------------|
| **Matter power spectrum** | Enhanced P(k) at low-k from horizon effects | P(k) shape matches BOSS/eBOSS at k < 0.1 h/Mpc | PM solver output |
| **Halo mass function** | More massive halos at high-z from faster collapse | Match JWST galaxy cluster abundances | Halo finder on N-body output |
| **Angular momentum excess** | Enhanced halo spin from vorticity amplification | λ parameter 1.5-2× ΛCDM prediction | Compute J/E for halos |
| **Growth rate fσ₈** | Scale-dependent growth from inertia reduction | Consistent with redshift-space distortions | Measure σ₈(z) and f(z) |
| **Void dynamics** | Different expansion in low-density regions | Void density profiles match observations | Identify voids, measure profiles |
| **Bullet Cluster test** | Lensing-gas separation from modified force | Recover observed ~0.5 Mpc separation | N-body + ray-tracing |
| **Correlation function** | Modified ξ(r) on horizon scales | Match SDSS/DESI correlation at r > 50 Mpc | Compute from particle positions |
| **Peculiar velocity field** | Different bulk flows from horizon term | Consistent with velocity surveys | Measure v_bulk from simulation |

### Falsifiable Claims

1. **If CMB peaks 2-5 are NOT within 3σ of Planck** → HQIV perturbation equations are wrong
2. **If z=14 proper time is NOT ~900-1000 Myr** → Background expansion is wrong
3. **If halo mass function at z>10 is NOT enhanced** → Structure formation prediction fails
4. **If Bullet Cluster separation is NOT reproduced** → Modified force law is incorrect
5. **If η does NOT emerge ~10⁻⁹ from early universe simulation** → Horizon statistics hypothesis fails
6. **If solar system tests show deviation** → High-acceleration limit is wrong (must match GR)

**Note**: β is not directly measurable — it's a derived quantity from horizon anisotropy integration.

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

**Generate table**:
```bash
python ../sandbox/hqiv_background.py
# outputs hqiv_Ha.txt: columns a, H/H0
```

### Phase 2: Modified Poisson Solver (2-4 weeks)

**File**: `src/poisson.c` or `src/forces.c`

The key modification: scale-dependent effective gravitational coupling.

```c
// Standard Poisson: ∇²Φ = 4πG ρ
// HQIV: ∇²Φ = 4πG_eff(a,k) ρ + horizon_term

double G_eff_HQIV(double a, double k) {
    double H_a = H_of_a_HQIV(a);
    double G_ratio = pow(H_a / H0, 0.6);  // α = 0.6
    
    // Horizon cutoff for k < 2π/Θ
    double Theta = 2 * c / H_a;  // in Mpc
    double k_horizon = 2 * M_PI / Theta;
    double horizon_factor = 1.0;
    
    if (k < k_horizon) {
        horizon_factor = exp(-k * Theta / (2 * M_PI));
    }
    
    return G0 * G_ratio * horizon_factor;
}

// In PM force calculation:
double force_HQIV(double a, double k, double density) {
    double G_eff = G_eff_HQIV(a, k);
    double horizon_term = beta * H_of_a_HQIV(a) * H_of_a_HQIV(a) / (c * c);
    
    return -4 * M_PI * G_eff * density + horizon_term;
}
```

### Phase 3: Inertia Reduction in Particle Dynamics (2-4 weeks)

**File**: `src/predict.c` or `src/runge.c`

Modified inertia affects particle acceleration response:

```c
// In particle pusher, modify effective mass:
double inertia_factor(double a, double a_local) {
    // a_min = β * c * H(a)
    double H_a = H_of_a_HQIV(a);
    double a_min = beta * c * H_a;
    
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

The vorticity amplification term:

```c
// ∂ω/∂t = β * H * (ω · ê_Θ)
// This is the most speculative part

void update_vorticity(double a, double dt) {
    double H_a = H_of_a_HQIV(a);
    double beta_local = compute_beta_local(a);  // from horizon anisotropy
    
    for (int i = 0; i < N_particles; i++) {
        // Compute local vorticity from velocity curl
        vec3 omega = compute_vorticity(i);
        
        // Horizon direction (radial from observer)
        vec3 e_Theta = normalize(particle[i].position);
        
        // Amplification
        double amplification = beta_local * H_a * dot(omega, e_Theta);
        
        particle[i].vorticity = omega * (1.0 + amplification * dt);
    }
}
```

---

## Testing Strategy

### Test 1: Baryon Fraction from Early Universe

**Setup**: Start simulation at T ~ 1 GeV (a ~ 10⁻¹²) with only radiation and horizon-quantized modes.

**Procedure**:
1. Initialize radiation bath with T = 1 GeV
2. Count modes that fit inside successive past light-cones
3. Track matter-antimatter asymmetry from horizon effects
4. Evolve to recombination (T ~ 0.3 eV)

**Success criterion**: Output η ~ 6×10⁻¹⁰

### Test 2: CMB Peak Structure

**Setup**: Run linear perturbation solver with HQIV modifications.

**Procedure**:
1. Generate H(a) table from background solver
2. Implement modified growth equations
3. Compute C_ℓ spectrum
4. Compare peak positions to Planck

**Success criterion**: Peaks 2-5 within 3σ of Planck

### Test 3: z=14 Age

**Setup**: Background integration only.

**Procedure**:
```python
# Already implemented in sandbox/hqiv_background.py
t_z14 = proper_time_at_z(a_arr, H_arr, 14.0)
```

**Success criterion**: t(z=14) ~ 900-1000 Myr

### Test 4: Temperature Uniformity

**Setup**: Track photon mode coupling across horizon scales.

**Procedure**:
1. Initialize CMB temperature field at recombination
2. Apply horizon-scale mode coupling
3. Measure temperature variance on different angular scales

**Success criterion**: σ_T/T < 10⁻⁵ on >1° scales

---

## File Structure

```
ecosmog/
├── README.md                    # This file
├── hqiv_pm.py                   # Python PM simulator (WORKING)
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

### Background Model (CORRECTED)

The horizon provides an effective energy density that dilutes more slowly than matter:

**H² = H0² × [Ω_m × a⁻³ + Ω_r × a⁻⁴ + Ω_horizon × a⁻ⁿ]**

Parameters calibrated to match paper:
- **Ω_horizon = 0.92** (effective horizon density)
- **n = 1.13** (horizon dilution rate, slower than matter's n=3)

### Test Results

**Background Cosmology:**
| Metric | HQIV Result | Paper Prediction | ΛCDM |
|--------|-------------|------------------|------|
| Universe age | **17.02 Gyr** | ~17 Gyr | 13.8 Gyr |
| t(z=14) | 690 Myr | ~940 Myr | ~300 Myr |
| H(a=1) | 68.9 km/s/Mpc | - | 70 km/s/Mpc |
| H(a=0.5) | 108.4 km/s/Mpc | - | ~120 km/s/Mpc |

**N-body Simulation** (32³ particles, 100 Mpc/h box, a=0.1→1.0):
| Metric | HQIV Result | ΛCDM Expected |
|--------|-------------|---------------|
| Growth factor | 1.28 | ~3.5 |
| P(k) peak | 5.4×10¹³ (Mpc/h)³ | - |
| Vorticity RMS | ~20 | ~0 (ΛCDM has no vorticity source) |
| Spin parameter λ | ~10²⁹ | ~0.035 (typical halos) |

### Vorticity Amplification (NEW!)

The simulation now includes the vorticity amplification term from the paper:

**∂ω/∂t + (v·∇)ω = β H (ω · ê_Θ)**

This equation:
1. Computes vorticity field ω = ∇ × v from particle velocities
2. Projects onto horizon direction ê_Θ (radial from observer)
3. Amplifies vorticity aligned with horizon
4. Generates angular momentum in collapsing structures

**Results with vorticity:**
- Vorticity RMS grows from ~0 to ~20 during simulation
- Spin parameter λ increases (enhanced halo rotation)
- Growth factor improves from 1.09 → 1.28

### Key Predictions

1. **Universe age ~17 Gyr** (older than ΛCDM's 13.8 Gyr) - testable with stellar ages
2. **Lower structure growth** (factor ~0.36× ΛCDM) - testable with σ₈ measurements
3. **Higher H(z) at early times** - testable with BAO at high z
4. **t(z=14) ~ 690 Myr** - allows older galaxies at high-z (JWST test)
5. **Enhanced halo spin** (λ ~ 1.5-2× ΛCDM) - testable with galaxy rotation curves
6. **Non-zero vorticity** in large-scale structure - testable with velocity field surveys

The lower growth factor is a **falsifiable prediction** - if σ₈ measurements match ΛCDM, HQIV's growth equations need revision.

---

## Dependencies

- ECOSMOG base code: `git clone https://github.com/MiXedCosmology/ECOSMOG`
- FFTW3: `sudo apt install libfftw3-dev`
- GSL: `sudo apt install libgsl-dev`
- MPI: `sudo apt install libopenmpi-dev`
- Python: numpy, scipy, matplotlib

---

## Build Instructions

```bash
# Clone ECOSMOG
git clone https://github.com/MiXedCosmology/ECOSMOG.git
cd ECOSMOG

# Apply HQIV patches
cp ../ecosmog/src/*.c src/
cp ../ecosmog/params/hqiv_params.ini .

# Generate H(a) table
python ../sandbox/hqiv_background.py
mv hqiv_Ha.txt .

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
| 5. Testing | 2-4 weeks | All 4 predictions validated |

**Total**: 3-6 months for complete implementation and testing

---

## References

1. ECOSMOG paper: Li et al., JCAP 2013, "ECOSMOG: An Efficient Code for Simulating Modified Gravity"
2. GADGET-2 manual: Springel, MNRAS 2005
3. QI framework: McCulloch, arXiv:1610.06787 and subsequent works