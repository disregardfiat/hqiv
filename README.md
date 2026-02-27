# Horizon-Quantized Informational Vacuum (HQIV)

Parameter-free cosmology from entanglement monogamy on causal horizons: zero dark energy, zero dark matter, horizon-driven late-time acceleration and first-principles baryon asymmetry \(\eta = 6.1\times10^{-10}\), spatial curvature \(\Omega_k^{\rm true}\approx +0.0098\), and fiducial global age **51.2 Gyr** (apparent ~13.8 Gyr from observer-centric lapse). Repo contains the paper, CLASS-HQIV integration, horizon lattice (baryogenesis + background table), and N-body Bullet test.

[Zenodo: Horizon-Quantized Informational Vacuum (HQIV): A Grand Unified Theory of Physics](https://zenodo.org/records/18717656)

---

## Quick start

### Paper (LaTeX)

From repo root:

```bash
cd paper
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Output: `paper/main.pdf`. Or upload `paper/` (main.tex + refs.bib + figures) to Overleaf.

### CLASS with HQIV (CMB / fiducial run)

The full CLASS tree is **not** in this repo (git-ignored). To reproduce the paper CMB figure and tables:

1. **Download CLASS**: `git clone https://github.com/lesgourg/class_public.git`
2. **Apply HQIV patches**: See **`class_hqiv_patches/README.md`** — copy sources/headers from `class_hqiv_patches/` into CLASS, add `hqiv.o` to the Makefile, then `make class`.
3. **Fiducial run**: Config `paper/class_fiducial_run.ini` (and mirror `class_hqiv_patches/paper_run/run.ini`). From the CLASS directory after building: `./class run.ini` (with `run.ini` copied from `class_hqiv_patches/paper_run/`).
4. **Peaks and CMB plot**: Run `class_hqiv_patches/extract_peaks.py` on the fiducial `*_cl.dat`; CMB figure: `class_hqiv_patches/plot_cmb_fiducial.py`. See paper § Run sets and reproducibility.

### Lattice table (baryogenesis → background for CLASS)

Produces the HQIV background table (emergent \(\Omega_m\), \(H_0\), \(\Omega_k^{\rm true}\), age) from the discrete light-cone:

```bash
cd horizon_modes/python
pip install numpy scipy   # optional: matplotlib for plots
python bulk.py            # or: from bulk import forward_4d_evolution; forward_4d_evolution(...)
```

Output: e.g. `hqiv_lattice_table.dat`. CLASS reads this with `hqiv_emergent = yes` and `hqiv_lattice_table = <path>`.

### Bullet Cluster N-body and lensing

Small test run (64³ particles, 20 steps, \(\gamma=0.40\)):

```bash
cd n-body_pysco_hqiv
pip install numpy scipy matplotlib   # + numba if available
python run_bullet.py --resolution 64 --npart 262144 --steps 20 --gamma 0.4 --box 5 --output ./output_bullet_64_20_g04
python postprocess_lensing.py --output ./output_bullet_64_20_g04 --box 5
```

Figure: `output_bullet_64_20_g04/lensing_comparison.png`. See **`n-body_pysco_hqiv/README.md`** for production runs and physics modules.

### HQVM algebra calculator (so(8) closure)

To verify full \(\mathfrak{so}(8)\) closure (28 dimensions) from \(\mathfrak{g}_2 + \Delta\) and confirm the hypercharge construction, use **`HQVM/matrices.py`** in a calculator or script:

```python
from HQVM.matrices import OctonionHQIVAlgebra
alg = OctonionHQIVAlgebra(verbose=False)
dim, history = alg.lie_closure_dimension()
if dim == 28:
    print("Full so(8) achieved — hypercharge construction valid")
```

Hypercharge verification (block entry error, eigenvalues of 4×4 block, and \([Y, \mathfrak{g}_2]\) check):

```python
alg = OctonionHQIVAlgebra(verbose=False)
c, Y, _ = alg.hypercharge_coefficients()
ver = alg.hypercharge_verify(Y)
print("Block entry error:", ver["block_entry_error"])
print("Eigenvalues imag:", ver["eigenvalues_i_block"])
print("Max [Y, g₂_colour] norm:", ver["max_commutation_with_g2"])   # ~1e-14 or better
```

Run from repo root (so that `HQVM` is on the Python path), or execute `python HQVM/matrices.py` from repo root for a quick closure check. Requires `numpy`. See **`HQVM/matrices.py`** for the full implementation (left-multiplication matrices \(L(e_i)\), phase-lift generator \(\Delta\), and Lie closure).

**Quantum Maxwell calculator (browser):** Open **`HQVM/quantum_maxwell_calculator.html`** in a browser for a self-contained app: phase-horizon Maxwell degrees of freedom (φ, ε(φ), μ(φ), δθ′), paper calculators (Higgs mass, curvature imprint δ_E(m), η, so(8) closure), interactive **Lie-closure visualiser** (iteration vs dimension 15→28, “Full so(8) achieved ✓” badge), **β-running engine** (sliders γ, T₀; α_EM, sin²θ_W, α_s down to M_Z; Table 1), Hypercharge Inspector, and multiple instances each plotting onto a shared graph.

**AI skill for HQIV physics:** The file **`skill-physics.md`** (at repo root) is a self-contained instruction set for AI agents. It summarizes the physics framework, links to the raw sources ([`HQVM/matrices.py`](https://raw.githubusercontent.com/disregardfiat/hqiv/main/HQVM/matrices.py), [`HQVM/quantum_maxwell_calculator.html`](https://raw.githubusercontent.com/disregardfiat/hqiv/main/HQVM/quantum_maxwell_calculator.html), [`paper/main.tex`](https://raw.githubusercontent.com/disregardfiat/hqiv/main/paper/main.tex)), and defines invariants and workflows that tools/agents should follow when modifying code or text.

---

## Repo layout (paths relative to repo root)

| Path | Description |
|------|-------------|
| **paper/** | LaTeX source: main.tex, refs.bib, class_fiducial_run.ini; figures (cmb_spectrum_fiducial.pdf, lensing_comparison.png, etc.) |
| **class_hqiv_patches/** | HQIV-modified CLASS files only; paper_run/run.ini; extract_peaks.py, plot_cmb_fiducial.py, background_cost_scan.py |
| **horizon_modes/python/** | Lattice baryogenesis and 4D table: bulk.py (forward_4d_evolution), discrete_baryogenesis_horizon.py |
| **n-body_pysco_hqiv/** | Bullet N-body: run_bullet.py, postprocess_lensing.py; hqiv_modifications/ (phi_field, g_eff, bullet_ic, etc.) |
| **HQVM/** | Octonion HQIV algebra: matrices.py (so(8) closure); quantum_maxwell_calculator.html (Phase-Horizon Maxwell + paper calculators, multi-instance graph) |
| **ecosmog/** | ECOSMOG / CLASS-HQIV run notes |

---

## Key numbers (from paper)

- **\(\eta = 6.1\times10^{-10}\)** (baryon asymmetry), **\(\Omega_k^{\rm true}\approx +0.0098\)** (true curvature), **\(\Omega_m = 0.0191\)** (fiducial baryon-only), **\(\gamma = 0.40\)** (thermodynamic overlap).
- **Global age**: 51.2 Gyr at fiducial cost minimum; apparent age ~13.8 Gyr from \(\phi\)-dependent lapse and time compression (~3.96×).
- **CMB**: Fiducial run gives peaks P1–P6 (Table in paper); peak-alignment cost ≈ 3.6; axis-of-evil band \(\ell\lesssim 20\).
- **Bullet**: 64³, 20 steps, \(\gamma=0.40\); lensing comparison in paper from `n-body_pysco_hqiv` run + postprocess_lensing.

---

## References

- Paper: **paper/main.tex** (and appendices: ADM lapse, variational horizon action, curvature-imprint normalization).
- CLASS build and run: **class_hqiv_patches/README.md**.
- N-body and lensing: **n-body_pysco_hqiv/README.md**.
