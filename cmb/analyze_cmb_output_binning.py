#!/usr/bin/env python3
"""Compare ordinary CMB bandpowers with angle-dependent multipole sums."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cmb.benchmark_cmb_ensemble import bin_spectrum, legendre_matrix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    folder = args.folder
    meta = json.loads((folder / "benchmark.json").read_text())
    if meta["nsims"] != 1:
        parser.error("this figure expects one paired realization")
    edges = np.asarray(meta["bin_edges"], dtype=int)
    lo, hi = int(edges[0]), int(edges[-1])
    fine_edges = np.append(np.arange(lo, hi, 5), hi)
    ell = np.arange(hi)

    with np.load(folder / "unbinned_spectra.npz") as data:
        key = "connected_q_residual"
        residual = data[f"pair_0_lensed_{key}"][:hi] - data[f"pair_0_unlensed_{key}"][:hi]

    angles = np.geomspace(2, 180, 450)
    example_angles = np.array([2, 10, 30])
    n_angles = len(angles)
    mu = np.cos(np.deg2rad(np.concatenate((angles, example_angles)) / 60))
    kernel = (2 * ell[:, None] + 1) * legendre_matrix(mu, hi - 1) / (4 * np.pi)
    per_mode = residual[lo:, None] * kernel[lo:]
    exact = per_mode.sum(axis=0)

    fig, axes = plt.subplots(3, 1, figsize=(9, 10))
    with (folder / "output_binning_check.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("scheme", "L_min", "L_max_inclusive", "L_eff_mode_weighted",
                         "paired_connected_q_bandpower", *(
                             f"exact_angular_contribution_at_{angle}_arcmin" for angle in example_angles
                         )))
        for bins, scheme, label, color in ((fine_edges, "width_5", r"$\Delta L=5$", "#6292b1"),
                                           (edges, "original", "Original bins", "#202c40")):
            centers = bin_spectrum(ell, bins)
            bandpowers = bin_spectrum(residual, bins)
            # Keep all Legendre kernels and approximate only the spectrum
            # within each bin. This isolates the effect of bin averaging.
            expanded = np.repeat(bandpowers, np.diff(bins))
            reconstructed = expanded @ kernel[lo:]
            error = reconstructed[:n_angles] - exact[:n_angles]
            axes[0].semilogx(centers, 1e7 * bandpowers, color=color, lw=1.1, label=label)
            axes[1].semilogx(angles, 1e4 * error, color=color, lw=1.2, label=label)
            print(f"{len(bins)-1} bands: largest angular error {np.max(np.abs(error)):.4g} "
                  f"({np.max(np.abs(error)) / np.max(np.abs(exact[:n_angles])):.2%} of the unbinned peak)")
            for lower, upper, center, value in zip(bins[:-1], bins[1:], centers, bandpowers):
                contribution = per_mode[lower-lo:upper-lo, n_angles:].sum(axis=0)
                writer.writerow((scheme, lower, upper - 1, center, value, *contribution))

    cumulative = np.cumsum(per_mode[:, n_angles:], axis=0)
    colors = ("#bb4b2f", "#9d659f", "#24718c")
    for angle, curve, color in zip(example_angles, cumulative.T, colors):
        axes[2].semilogx(ell[lo:], 1e4 * curve, color=color, lw=1.35, label=fr"${angle}'$")
    for limit in (100, 600):
        if lo < limit < hi:
            axes[2].axvline(limit, color="0.55", lw=0.8, ls="--")
    axes[0].set(title="Paired connected residual in multipole bins",
                xlabel=r"Mode-weighted bin center $L_{\rm eff}$", ylabel=r"$10^7\,\Delta T_{q,L}$")
    axes[1].set(title="Error from replacing the spectrum within each bin by its mean",
                xlabel=r"Separation $\theta$ [arcmin]", ylabel=r"$10^4\,[\Delta T_{q,\rm binned}-\Delta T_{q,\rm exact}]$")
    axes[2].set(title="Angular residual accumulated from individual multipoles",
                xlabel=r"Highest included output multipole $L$", ylabel=r"$10^4\,\Delta T_q(\theta)$")
    axes[0].set_xlim(lo, min(hi, 1000))
    axes[1].set_xlim(angles[0], angles[-1])
    axes[2].set_xlim(lo, hi - 1)
    for ax in axes:
        ax.axhline(0, color="0.55", lw=0.8)
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, ncol=2)
    fig.suptitle(fr"CMB lensing pair with input temperature $\ell\geq{meta['temperature_lmin']}$", fontsize=13)
    fig.tight_layout()
    fig.savefig(folder / "output_binning_check.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    main()
