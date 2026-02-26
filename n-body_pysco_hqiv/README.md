# HQIV N-body Simulation with PySCo

This module implements the full action-derived HQIV (Horizon-Quantized Informational Vacuum) physics for N-body simulations, specifically targeting the Bullet Cluster test case.

## Paper Reference

All equations are derived from:
- **paper/main.tex**: Main HQIV framework paper
- **paper/maxwell.tex**: Derivation from Maxwell's equations and Schuller's geometric framework

## Key Physics

### Modified Einstein Equation (paper/main.tex, Section 4)
```
G_μν + γ(φ/c²) g_μν = (8π G_eff(φ)/c⁴) T_μν
```

### Background Friedmann Equation (Section 5)
```
3H² - γH = 8π G_eff(H) (ρ_m + ρ_r)
```

### Inertia Factor (Section 4.2)
```
f(a_loc, φ) = max( a_loc / (a_loc + cφ/6), f_min )
```

### Vorticity Source (Section 6)
```
∂ω/∂t + (v·∇)ω = (∂f/∂φ) (k × ∇φ) · ê_ω
```

### Varying Gravitational Coupling (Section 3)
```
G(a) = G0 × (H(a)/H0)^α
```

## Directory Structure

```
n-body_pysco_hqiv/
├── README.md                    # This file
├── config_hqiv.ini             # Configuration (copied to hqiv_modifications/)
├── run_bullet.py               # Main simulation script
├── postprocess_lensing.py      # Weak lensing post-processing
├── hqiv_modifications/         # HQIV physics modules
│   ├── __init__.py
│   ├── config_hqiv.ini        # Configuration file
│   ├── phi_field.py           # Geometric horizon field φ(x)
│   ├── inertia_factor.py      # Thermodynamic inertia reduction
│   ├── vorticity_source.py    # Action-derived vorticity source
│   ├── g_eff.py               # Varying gravitational coupling
│   └── bullet_ic.py           # Bullet Cluster initial conditions
├── pysco/                      # PySCo submodule (git submodule)
└── tests/
    └── test_hqiv.py           # Unit tests
```

## Installation

```bash
# Clone the repository with submodules
git clone --recursive https://github.com/your-repo/HQIV.git
cd HQIV/n-body_pysco_hqiv

# Install dependencies
pip install numpy scipy matplotlib

# Install PySCo (if not using submodule)
pip install pysco-cosmo
```

## Quick Start

### Run a small test simulation

```bash
# Run with default parameters (small test)
python run_bullet.py --resolution 64 --npart 1e5 --steps 100

# Run with custom configuration
python run_bullet.py --config hqiv_modifications/config_hqiv.ini --resolution 128 --npart 1e6
```

### Long / production runs (detached + resource limits)

Heavy runs (e.g. 512³, 10⁷ particles) can **exhaust RAM/CPU and kill Cursor** if run on the same machine without limits. Prefer one of:

1. **Run on another machine** or a job queue (SLURM, etc.), or  
2. **Cap resources** so the run cannot starve the rest of the system:

```bash
# With memory and CPU caps (recommended when running on same box as Cursor)
BULLET_MEM_GB=8 BULLET_CPU_PCT=70 ./run_bullet_nohup.sh ./bullet_5mpc_512/ bullet_5mpc_512.log
```

Without `BULLET_MEM_GB` / `BULLET_CPU_PCT`, the script still runs detached but prints a warning; the run can freeze or kill the IDE.

```bash
# Custom output dir and log
./run_bullet_nohup.sh ./my_run/ my_run.log

# Extra run_bullet.py args (after output dir and log)
./run_bullet_nohup.sh ./my_run/ my_run.log --steps 300 --gamma 0.40
```

Manual detached run: `nohup python3 -u run_bullet.py ... >> run.log 2>&1 & disown -h $!`

### Post-process results for lensing

```bash
# Compute convergence map from final snapshot
python postprocess_lensing.py --output ./output/ --snapshot final_particles.npz

# Use full geodesic ray-tracing
python postprocess_lensing.py --output ./output/ --method full
```

### Run tests

```bash
# Run all tests
python tests/test_hqiv.py

# Or with pytest
python -m pytest tests/test_hqiv.py -v
```

## Configuration

Edit `hqiv_modifications/config_hqiv.ini` to modify parameters:

### Cosmology Parameters
```ini
[cosmology]
H0 = 73.2              # Hubble constant [km/s/Mpc]
Omega_m = 0.06         # Matter density (baryons only in HQIV)
Neff = 3.046           # Effective neutrino species
T_CMB = 2.7255         # CMB temperature [K]
```

### HQIV Parameters
```ini
[hqiv_parameters]
gamma = 0.40           # Thermodynamic coefficient (Brodie's overlap integral)
alpha_G = 0.6          # Exponent for varying G
chi = 0.172            # Light-cone average scaling factor
f_min = 0.01           # Thermodynamic floor for inertia
inertia_form = thermo  # 'thermo', 'sqrt', or 'brodie'
```

### Bullet Cluster Parameters
```ini
[bullet_cluster]
z_obs = 0.3            # Observation redshift
M_main = 2.5e14        # Main cluster mass [Msun]
M_sub = 1.5e14         # Subcluster mass [Msun]
v_collision = 4500     # Collision velocity [km/s]
impact_param = 150     # Impact parameter [kpc]
gas_fraction = 0.15    # Gas fraction
```

## Output

The simulation produces:

1. **Snapshots**: `snapshot_XXXXX.npz` - Particle positions, velocities, masses at each save interval
2. **Final results**: `final_particles.npz` - Final particle data
3. **History**: `simulation_history.npz` - Time evolution of diagnostics
4. **Density maps**: `density_maps.npz` - 3D density field

Post-processing produces:

1. **Convergence map**: κ field for weak lensing
2. **Shear field**: γ₁, γ₂ components
3. **Magnification map**: μ field
4. **Comparison plot**: `lensing_comparison.png`

## Key Predictions

The HQIV framework makes several testable predictions for the Bullet Cluster:

1. **Gas-Galaxy Offset**: ~180 kpc separation between X-ray gas and lensing peak
2. **Lensing/Gas Ratio**: ~4-6:1 ratio of lensing mass to gas mass
3. **Vorticity Growth**: Positive growth exponent (β ≈ +1.9) vs. standard decay (β = -2)
4. **G_eff Evolution**: G_eff > G0 at high redshift, affecting early structure formation

## Physics Modules

### phi_field.py
Computes the geometric horizon field φ(x) = 2c²/Θ_local(x) from:
- Expansion scalar θ = ∇·v
- Local horizon radius Θ_local
- Background H(a) from modified Friedmann equation

### inertia_factor.py
Implements the thermodynamic inertia reduction:
- f(a_loc, φ) = max(a_loc / (a_loc + cφ/6), f_min)
- Species-dependent factors (gas vs. galaxies)
- Direction-dependent effects

### vorticity_source.py
Computes the action-derived vorticity source:
- Source term: (∂f/∂φ) (k × ∇φ) · ê_ω
- Vorticity evolution with amplification
- Growth exponent calculation

### g_eff.py
Implements varying gravitational coupling:
- G_eff(a) = G0 (H(a)/H0)^α
- Horizon correction to Poisson equation
- Modified force calculation

### bullet_ic.py
Generates Bullet Cluster initial conditions:
- Two NFW halos with specified masses and concentrations
- Gas particles with thermal distribution
- Galaxy particles following dark matter profile

## Comparison with Observations

The Bullet Cluster observations (Clowe et al. 2006) show:
- Offset between X-ray gas and weak lensing peak: ~180 kpc
- Lensing/gas mass ratio: ~4-6:1
- Collision velocity: ~4500 km/s

HQIV should reproduce these without dark matter through:
- Direction-dependent inertia reduction
- Different effective masses for gas vs. galaxies
- Modified gravitational dynamics

## References

1. Ettinger, S. Jr. (2026). "Horizon-Quantized Informational Vacuum (HQIV): A Covariant Baryon-Only Cosmological Framework from Quantised Inertia"
2. Brodie (2026). Thermodynamic derivation of MOND-like modification
3. Clowe et al. (2006). "A Direct Empirical Proof of the Existence of Dark Matter"
4. McCulloch, M.E. (2007-2016). Quantised Inertia papers

## License

MIT License - See LICENSE file for details.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Contact

For questions or issues, please open a GitHub issue or contact the HQIV Team.