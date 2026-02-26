# Bullet 5 Mpc, 512³ Lensing Report

- Output directory: `/home/jr/Repos/HQIV/n-body_pysco_hqiv/bullet_5mpc_512`
- Ray-tracing method: `born`
- Snapshots: snapshot_00010.npz, snapshot_00020.npz, snapshot_00030.npz, snapshot_00040.npz

## Summary Table

| Snapshot | a | z | κ_min | κ_max | n_peaks | offset (lensing–gas) [kpc] |
|----------|---|---|-------|-------|---------|-----------------------------|
| snapshot_00010.npz | 0.510000 | 0.961 | 0.0000e+00 | 8.8835e-02 | 4 | 2505.2 |
| snapshot_00020.npz | 0.520000 | 0.923 | 0.0000e+00 | 9.1237e-02 | 2 | 568.3 |
| snapshot_00030.npz | 0.530000 | 0.887 | 0.0000e+00 | 9.1740e-02 | 3 | 221.1 |
| snapshot_00040.npz | 0.540000 | 0.852 | 0.0000e+00 | 8.9349e-02 | 2 | 568.7 |

## Figures

- `kappa_evolution_4snapshots.png`: κ evolution across the four snapshots.
- `offset_vs_redshift.png`: Evolution of the lensing–gas offset as a function of redshift.
- `lensing_snapshot_*.png`: Per-snapshot κ, gas projection, magnification, and shear.

All figures live in the same directory as this report and can be directly included in LaTeX or other writeups.