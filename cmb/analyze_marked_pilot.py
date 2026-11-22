#!/usr/bin/env python3
"""Show the marked-spectrum response across paired CMB seeds."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cmb.analyze_cmb_high_ell import correlation
from cmb.benchmark_cmb_ensemble import bin_spectrum


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    folder = args.folder
    meta = json.loads((folder / "benchmark.json").read_text())
    edges = np.asarray(meta["bin_edges"], dtype=int)
    centers = bin_spectrum(np.arange(meta["lmax"] + 1), edges)
    theta = np.geomspace(2, 180, 320)

    paired_cl = []
    lensed_cl = []
    unlensed_cl = []
    paired_xi = []
    with np.load(folder / "unbinned_spectra.npz") as saved:
        for seed_index in range(meta["nsims"]):
            lensed = saved[f"pair_{seed_index}_lensed_marked_residual"]
            unlensed = saved[f"pair_{seed_index}_unlensed_marked_residual"]
            paired = lensed - unlensed
            lensed_cl.append(bin_spectrum(lensed, edges))
            unlensed_cl.append(bin_spectrum(unlensed, edges))
            paired_cl.append(bin_spectrum(paired, edges))
            paired_xi.append(correlation(paired, theta, int(edges[0]), int(edges[-1])))

    lensed_cl = np.asarray(lensed_cl)
    unlensed_cl = np.asarray(unlensed_cl)
    paired_cl = np.asarray(paired_cl)
    paired_xi = np.asarray(paired_xi)
    fig, axes = plt.subplots(3, 1, figsize=(9, 10))
    colors = plt.cm.viridis(np.linspace(0.18, 0.78, meta["nsims"]))
    for seed_index, color in enumerate(colors):
        axes[0].plot(centers, paired_cl[seed_index] * 1e7, color=color, lw=1,
                     label=f"seed {seed_index + 1}")
        axes[2].semilogx(theta, paired_xi[seed_index] * 1e4, color=color, lw=1)
    axes[0].plot(centers, paired_cl.mean(axis=0) * 1e7, color="#202c40",
                 lw=2.4, label="mean")
    axes[1].plot(centers, lensed_cl.mean(axis=0) * 1e7, color="#bb4b2f",
                 lw=1.7, label="lensed mean")
    axes[1].plot(centers, unlensed_cl.mean(axis=0) * 1e7, color="#24718c",
                 lw=1.7, label="unlensed mean")
    axes[2].semilogx(theta, paired_xi.mean(axis=0) * 1e4, color="#202c40",
                     lw=2.4, label="mean")
    axes[0].set(title="Marked power spectrum, lensed minus unlensed after phase controls",
                xlabel=r"Output multipole $L$", ylabel=r"$10^7\Delta C_L^{\rm mark}$")
    axes[1].set(title="Residual marked spectra for each map class",
                xlabel=r"Output multipole $L$", ylabel=r"$10^7 C_L^{\rm mark,res}$")
    axes[2].set(title="Paired marked angular correlation from the same output range",
                xlabel=r"Separation $\theta$ [arcmin]", ylabel=r"$10^4\Delta \xi_m(\theta)$")
    for ax in axes:
        ax.axhline(0, color="0.45", lw=0.8)
        ax.grid(alpha=0.16)
        ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.suptitle(f"{meta['nsims']} paired skies, {meta['surrogates']} phase controls per map",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(folder / "marked_pilot.png", dpi=190)
    plt.close(fig)

    with (folder / "marked_pilot_summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("seed_index", "kind", "L_min", "L_max", "marked_bandpower"))
        for seed_index in range(meta["nsims"]):
            for kind, values in (("lensed", lensed_cl[seed_index]),
                                 ("unlensed", unlensed_cl[seed_index]),
                                 ("paired", paired_cl[seed_index])):
                for lower, upper, value in zip(edges[:-1], edges[1:], values):
                    writer.writerow((seed_index, kind, lower, upper - 1, value))

    low, high = 100, min(600, int(edges[-1]))
    if high > low:
        with np.load(folder / "unbinned_spectra.npz") as saved:
            for seed_index in range(meta["nsims"]):
                paired = (saved[f"pair_{seed_index}_lensed_marked_residual"]
                          - saved[f"pair_{seed_index}_unlensed_marked_residual"])
                band = bin_spectrum(paired, np.array([low, high]))[0]
                short_angle = correlation(paired, np.array([2.0]), low, high)[0]
                print(f"seed {seed_index + 1}: paired marked C_L[{low},{high})={band:.3e}; "
                      f"xi_m(2 arcmin)={short_angle:.3e}")


if __name__ == "__main__":
    main()
