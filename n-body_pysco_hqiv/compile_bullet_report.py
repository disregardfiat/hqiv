#!/usr/bin/env python3
"""
Compile Bullet 5 Mpc, 512³ run snapshots into lensing graphics and a short report.

This script processes the four stored snapshots in `bullet_5mpc_512/` using the
existing `postprocess_lensing.LensingPostprocessor` and generates:

1. Per-snapshot comparison plots:
   - `lensing_snapshot_00010.png`, ..., `lensing_snapshot_00040.png`
2. A 2×2 panel showing κ for all four snapshots:
   - `kappa_evolution_4snapshots.png`
3. An offset vs redshift plot:
   - `offset_vs_redshift.png`
4. A markdown report summarizing key diagnostics:
   - `bullet_5mpc_512_report.md`
"""

from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np

from postprocess_lensing import LensingPostprocessor


def process_snapshots(
    output_dir: Path,
    snapshots: List[str],
    method: str = "born",
) -> List[dict]:
    """Run lensing post-processing on a list of snapshots."""
    processor = LensingPostprocessor(output_dir)
    results_list: List[dict] = []

    for snap in snapshots:
        results = processor.process(snap, method=method)

        # Per-snapshot comparison graphic
        png_name = f"lensing_{snap.replace('.npz', '')}.png"
        processor.create_comparison_plot(results, output_file=png_name)

        results_list.append(results)

    return results_list


def make_kappa_evolution_figure(
    output_dir: Path,
    results_list: List[dict],
    snapshots: List[str],
) -> None:
    """Create a 2×2 panel of κ maps for the four snapshots."""
    if not results_list:
        return

    box_size = float(results_list[0]["kappa"].shape[0])
    extent = [0, box_size, 0, box_size]

    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    axes = axes.flatten()

    for ax, res, snap in zip(axes, results_list, snapshots):
        kappa = res["kappa"]
        z = 1.0 / res["a"] - 1.0
        im = ax.imshow(
            kappa.T,
            origin="lower",
            extent=extent,
            cmap="viridis",
        )
        ax.set_title(f"{snap} (z={z:.3f})")
        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")
        fig.colorbar(im, ax=ax, label="κ")

    plt.tight_layout()
    out_path = output_dir / "kappa_evolution_4snapshots.png"
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


def make_offset_vs_redshift_figure(
    output_dir: Path,
    results_list: List[dict],
    snapshots: List[str],
) -> None:
    """Plot lensing–gas offset as a function of redshift."""
    if not results_list:
        return

    zs = [1.0 / res["a"] - 1.0 for res in results_list]
    offsets = [res["offset_kpc"] for res in results_list]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(zs, offsets, "o-", label="Offset (lensing–gas)")

    for z, off, snap in zip(zs, offsets, snapshots):
        ax.annotate(
            snap.replace(".npz", ""),
            (z, off),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=8,
        )

    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Offset [kpc]")
    ax.set_title("Bullet Cluster offset vs redshift (HQIV)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    out_path = output_dir / "offset_vs_redshift.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


def write_markdown_report(
    output_dir: Path,
    results_list: List[dict],
    snapshots: List[str],
    method: str,
) -> None:
    """Write a compact markdown report summarizing the four snapshots."""
    if not results_list:
        return

    lines: List[str] = []
    lines.append("# Bullet 5 Mpc, 512³ Lensing Report")
    lines.append("")
    lines.append(f"- Output directory: `{output_dir}`")
    lines.append(f"- Ray-tracing method: `{method}`")
    lines.append(f"- Snapshots: {', '.join(snapshots)}")
    lines.append("")
    lines.append("## Summary Table")
    lines.append("")
    lines.append(
        "| Snapshot | a | z | κ_min | κ_max | n_peaks | offset (lensing–gas) [kpc] |"
    )
    lines.append(
        "|----------|---|---|-------|-------|---------|-----------------------------|"
    )

    for snap, res in zip(snapshots, results_list):
        a = float(res["a"])
        z = 1.0 / a - 1.0
        kmin = float(res["kappa"].min())
        kmax = float(res["kappa"].max())
        n_peaks = len(res["peaks"])
        offset = float(res["offset_kpc"])
        lines.append(
            f"| {snap} | {a:.6f} | {z:.3f} | {kmin:.4e} | {kmax:.4e} | {n_peaks:d} | {offset:.1f} |"
        )

    lines.append("")
    lines.append("## Figures")
    lines.append("")
    lines.append("- `kappa_evolution_4snapshots.png`: κ evolution across the four snapshots.")
    lines.append(
        "- `offset_vs_redshift.png`: Evolution of the lensing–gas offset as a function of redshift."
    )
    lines.append(
        "- `lensing_snapshot_*.png`: Per-snapshot κ, gas projection, magnification, and shear."
    )
    lines.append("")
    lines.append(
        "All figures live in the same directory as this report "
        "and can be directly included in LaTeX or other writeups."
    )

    report_path = output_dir / "bullet_5mpc_512_report.md"
    report_path.write_text("\n".join(lines))


def main() -> None:
    base_dir = Path(__file__).parent
    output_dir = base_dir / "bullet_5mpc_512"

    snapshots = [
        "snapshot_00010.npz",
        "snapshot_00020.npz",
        "snapshot_00030.npz",
        "snapshot_00040.npz",
    ]
    method = "born"

    results_list = process_snapshots(output_dir, snapshots, method=method)

    make_kappa_evolution_figure(output_dir, results_list, snapshots)
    make_offset_vs_redshift_figure(output_dir, results_list, snapshots)
    write_markdown_report(output_dir, results_list, snapshots, method=method)


if __name__ == "__main__":
    main()

