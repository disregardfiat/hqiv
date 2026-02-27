#!/usr/bin/env bash
# One-click reproducibility for paper figures and tables (see main.tex § Run sets and reproducibility).
# Run from repository root: ./paper/reproduce_figures.sh
# Optional: set CLASS_DIR to the CLASS build directory (class_public) to run CLASS + CMB figure.

set -e
HQIV_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HQIV_ROOT"

echo "=== 1. Lattice + beta engine (precision constants) ==="
cd "$HQIV_ROOT/horizon_modes/python"
python3 beta_engine.py --n_steps 6000 --out "$HQIV_ROOT/paper/precision_constants.dat"
cd "$HQIV_ROOT"
echo "  -> paper/precision_constants.dat, horizon_modes/python/hqiv_lattice_table.dat"

if [ -n "${CLASS_DIR}" ] && [ -f "${CLASS_DIR}/class" ]; then
  echo "=== 2. CLASS fiducial run ==="
  cp -f "$HQIV_ROOT/class_hqiv_patches/paper_run/run.ini" "${CLASS_DIR}/"
  ( cd "${CLASS_DIR}" && ./class run.ini )
  # Copy outputs so extract_peaks and plot_cmb_fiducial find them (root = hqiv_min)
  cp -f "${CLASS_DIR}"/hqiv_min*.dat "$HQIV_ROOT/class_hqiv_patches/paper_run/" 2>/dev/null || true
  echo "=== 3. Extract peaks + CMB figure ==="
  python3 "$HQIV_ROOT/class_hqiv_patches/extract_peaks.py" "$HQIV_ROOT/class_hqiv_patches/paper_run/hqiv_min00_cl.dat"
  python3 "$HQIV_ROOT/class_hqiv_patches/plot_cmb_fiducial.py"
  echo "  -> paper/cmb_spectrum_fiducial.pdf, peaks in class_hqiv_patches/paper_run/"
else
  echo "=== 2–3. CLASS + CMB figure skipped (set CLASS_DIR to CLASS build dir to run)"
  echo "  Example: CLASS_DIR=/path/to/class_public ./paper/reproduce_figures.sh"
fi

echo "Done."
