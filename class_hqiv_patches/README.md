# CLASS-HQIV patches

This directory contains **only the HQIV-modified files** for the [CLASS](http://class-code.net) Boltzmann code. The full CLASS tree and the full `class_hqiv` fork are **not** tracked in this repository (they are in `.gitignore`). To build CLASS with HQIV support you must **download CLASS and apply these patches**.

## 1. Download CLASS

Clone the official CLASS repository:

```bash
git clone https://github.com/lesgourg/class_public.git
cd class_public
```

(Or download a release tarball from the [CLASS website](http://class-code.net).)

## 2. Apply HQIV patches

From the **CLASS root** (the `class_public` directory):

```bash
# Overwrite/copy HQIV source and headers
cp /path/to/HQIV/class_hqiv_patches/source/hqiv.c       source/
cp /path/to/HQIV/class_hqiv_patches/source/background.c source/
cp /path/to/HQIV/class_hqiv_patches/source/input.c      source/
cp /path/to/HQIV/class_hqiv_patches/source/perturbations.c source/
cp /path/to/HQIV/class_hqiv_patches/source/thermodynamics.c source/
cp /path/to/HQIV/class_hqiv_patches/include/hqiv.h      include/
cp /path/to/HQIV/class_hqiv_patches/include/background.h include/
```

Replace `/path/to/HQIV` with the path to this HQIV repository (or use relative paths if you cloned HQIV next to `class_public`).

### Makefile change

Add `hqiv.o` to the `SOURCE` line in the CLASS `Makefile`. Find the line that looks like:

```makefile
SOURCE = input.o background.o thermodynamics.o perturbations.opp primordial.opp fourier.o transfer.opp harmonic.opp lensing.opp distortions.o
```

and change it to:

```makefile
SOURCE = input.o background.o thermodynamics.o perturbations.opp primordial.opp fourier.o transfer.opp harmonic.opp lensing.opp distortions.o hqiv.o
```

## 3. Build and run

```bash
make clean
make class
./class explanatory.ini   # test vanilla run first
```

Copy the HQIV test parameter files from `class_hqiv_patches/` (e.g. `test_hqiv.ini`, `test_minimal.ini`) into the CLASS directory, then:

```bash
./class test_hqiv.ini
```

Output files use the prefix set by `root` in the `.ini` file (e.g. `hqiv_background.dat`).

### 4. Peak-alignment input-space search

With the **Omega_eff closure** fix, the budget equation sets `Omega_Lambda = 1 - γ/3 - Ω_tot` so that `H(a=1) = H0`. Universe age then matches the Python wrapper (~19 Gyr for γ=0.4, not ~34 Gyr).

To re-run the multipole peak alignment search (paper: *peak_alignment_scan*):

```bash
# From CLASS root (class_public), with Python interface built:
make
export PYTHONPATH=$PWD/python:$PYTHONPATH
python /path/to/HQIV/class_hqiv_patches/peak_alignment_scan.py --steps 4 --out peak_scan_results
```

Optional: `--fine` runs a second pass around the best point. Results are written to `peak_scan_results.txt` and `peak_scan_results.npy`.

## Contents of this directory

| Path | Description |
|------|-------------|
| `include/hqiv.h` | HQIV API and parameters |
| `include/background.h` | Background module header (HQIV-modified) |
| `source/hqiv.c` | HQIV background and \(G_{\rm eff}\) logic |
| `source/background.c` | Modified Friedmann \(3H^2 - \gamma H = 8\pi G_{\rm eff}\,\rho\) |
| `source/input.c` | Parser for `HQIV`, `HQIV_gamma`, `HQIV_alpha`, etc. |
| `source/perturbations.c` | Perturbation equations with HQIV terms |
| `source/thermodynamics.c` | Thermodynamics with HQIV hooks |
| `test_hqiv.ini` | Example HQIV run (baryon-only, \(\gamma=0.4\)) |
| `test_minimal.ini` | Minimal HQIV parameters |
| `test_hqiv_debug.ini` | Debug/verbose run |
| `test_baryons_only.ini` | Baryon-only test |
| `peak_alignment_scan.py` | Input-space search over γ, ω_b, h, α to match CMB peak positions to Planck |
| `paper_run/run.ini` | Fiducial run config (mirror of `paper/class_fiducial_run.ini`); run from CLASS build dir to reproduce paper CMB/tables |
| `extract_peaks.py` | Extract P1–P6 peak positions from fiducial `*_cl.dat` |
| `plot_cmb_fiducial.py` | Produce CMB TT plot (paper figure) |
| `background_cost_scan.py` | Peak-alignment cost at fiducial point |

For more run options and parameter descriptions, see **Running the CLASS-HQIV code** in `ecosmog/README.md`.
