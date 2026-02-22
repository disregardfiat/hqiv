# Horizon-Quantized Informational Vacuum (HQIV)

Parameter-minimal cosmology derived from **Quantised Inertia** (McCulloch 2007–2026): zero dark energy, zero dark matter, horizon-driven late-time acceleration and faster early structure formation. Repo contains the arXiv-ready paper, sandbox numerics, and CLASS/HiCLASS integration scaffolding.

[Zenodo Preprint: Horizon-Quantized Informational Vacuum (HQIV): A Covariant Baryon-Only Cosmological Framework from Quantised Inertia](https://zenodo.org/records/18717656)

---

## Quick start

### Paper (LaTeX)

```bash
cd paper
pdflatex main
bibtex main
pdflatex main
pdflatex main
# Or use Overleaf: upload paper/ (main.tex + refs.bib).
```

### Sandbox (background integration)

Produces universe age, proper time at z=14, d_A(z=14), and `hqiv_Ha.txt` for CLASS.

**Recommended: fresh venv with scipy** (needed for stable ODE integration; numpy-only RK4 can be unstable):

```bash
python -m venv hqiv_env
source hqiv_env/bin/activate   # or hqiv_env\Scripts\activate on Windows
pip install numpy scipy matplotlib classy
cd sandbox
python hqiv_background.py
```

If you only have numpy, the script falls back to a built-in RK4; the run may then give off numbers until scipy is available. Paper targets: **17.2 Gyr** age, **940 Myr** at z=14, **d_A(z=14) ≈ 1.42 Gpc**.

**CLASS from source (for custom background / HiCLASS-style modified gravity):**

```bash
git clone https://github.com/lesgourg/class_public.git
cd class_public
make -j4
pip install -e .
```

(For the modified-gravity fork with μ/γ support, use the HiCLASS branch; see `class_integration/MODIFICATIONS.md`.)

### CLASS with HQIV (Boltzmann / CMB)

The full CLASS tree and the CLASS-HQIV fork are **not** tracked in this repo (they are git-ignored). To build and run CLASS with HQIV support:

1. **Download CLASS**: `git clone https://github.com/lesgourg/class_public.git`
2. **Apply HQIV patches**: Copy the files from `class_hqiv_patches/` into the CLASS `source/` and `include/` directories, and add `hqiv.o` to the Makefile `SOURCE` line. Full steps are in **`class_hqiv_patches/README.md`**.
3. **Build and run**: `make class` then `./class test_hqiv.ini` (copy the `.ini` files from `class_hqiv_patches/`).

See **`ecosmog/README.md`** (section “Running the CLASS-HQIV code”) for parameter descriptions and run options.

---

## Repo layout

| Path | Contents |
|------|----------|
| `paper/` | LaTeX source (main.tex, refs.bib) — arXiv/Overleaf ready |
| `sandbox/` | HQIV background ODE (scipy), H(a) table output |
| `class_hqiv_patches/` | **HQIV-modified CLASS files only** — clone CLASS and apply these to build CLASS-HQIV |
| `ecosmog/` | ECOSMOG implementation plan and **steps to run CLASS-HQIV** |
| `class_integration/` | CLASS/HiCLASS params, H(a) export (legacy) |

---

## Key numbers (from paper)

- **β ≈ 0.81**, visible **Ω_m = 0.048** (two effective parameters, to be eliminated).
- Age today: **17.2 Gyr**; proper time at z=14: **940 Myr** (≈3.2× ΛCDM).
- CMB: first peak ℓ ≈ 217; low-ℓ power suppressed ≈18%; χ²_SN+BAO better than untuned ΛCDM by Δχ² = −14.5 (2 fewer parameters).
- Redshift quantization clusters at z ≈ 10.2, 11.8, 13.5, 14.9 (JWST-testable).


