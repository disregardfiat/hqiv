#!/usr/bin/env python3
"""
Extract CMB TT acoustic peak positions (P1--P6) from a CLASS *_cl.dat file.
D_ell = l(l+1)C_l/(2pi) is in column 2. Use window-based search so we always
report 6 peaks (max in each of 6 multipole windows).
"""
import os
import sys
import numpy as np

# Approximate Planck peak positions; we search in windows around these
PEAK_WINDOWS = [(50, 400), (400, 650), (650, 950), (950, 1350), (1350, 1550), (1550, 2500)]

def main():
    cl_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(__file__), 'paper_run', 'hqiv_min00_cl.dat')
    if not os.path.isfile(cl_file):
        print("Usage: extract_peaks.py [path/to/*_cl.dat]", file=sys.stderr)
        sys.exit(1)
    data = np.loadtxt(cl_file, comments='#')
    if data.ndim == 1:
        data = data.reshape(-1, 2)
    ell = data[:, 0]
    d_ell = data[:, 1]
    peaks_ell = []
    for lo, hi in PEAK_WINDOWS:
        mask = (ell >= lo) & (ell <= hi)
        if not np.any(mask):
            continue
        idx = np.argmax(d_ell[mask])
        ell_win = ell[mask]
        peaks_ell.append(int(round(ell_win[idx])))
    for j, L in enumerate(peaks_ell, 1):
        print("P{} = {}".format(j, L))
    return 0

if __name__ == '__main__':
    sys.exit(main())
