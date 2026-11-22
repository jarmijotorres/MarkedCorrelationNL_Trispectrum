#!/usr/bin/env python3
"""Compare paired output spectra for two input-temperature high-pass cuts."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cmb.benchmark_cmb_ensemble import bin_spectrum


FOLDERS = {
    100: Path("archive/cmb_pilots/run_cmb_actlike_lmin100_2048"),
    600: Path("archive/cmb_pilots/run_cmb_actlike_high_ell_2048"),
}
OUT = FOLDERS[100]


def main() -> None:
    results = {}
    for cut, folder in FOLDERS.items():
        meta = json.loads((folder / "benchmark.json").read_text())
        if meta["temperature_lmin"] != cut:
            raise ValueError(f"Unexpected temperature cut in {folder}")
        edges = np.asarray(meta["bin_edges"], dtype=int)
        with np.load(folder / "unbinned_spectra.npz") as saved:
            paired = {
                name: saved[f"pair_0_lensed_{name}"] - saved[f"pair_0_unlensed_{name}"]
                for name in ("q_phi", "marked_residual", "q_component_residual")
            }
        results[cut] = {name: bin_spectrum(values, edges) for name, values in paired.items()}
    centers = (edges[:-1] + edges[1:] - 1) / 2
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for cut, color in ((100, "#bb4b2f"), (600, "#24718c")):
        axes[0].plot(centers, 1e7 * centers * (centers + 1) * results[cut]["q_phi"],
                     color=color, lw=1.4, marker=".", ms=3, label=fr"$\ell_{{\min}}={cut}$")
        axes[1].plot(centers, 1e7 * results[cut]["marked_residual"],
                     color=color, lw=1.4, marker=".", ms=3, label=fr"$\ell_{{\min}}={cut}$")
        axes[2].plot(centers, 1e7 * results[cut]["q_component_residual"],
                     color=color, lw=1.4, marker=".", ms=3, label=fr"$\ell_{{\min}}={cut}$")
    axes[0].set_title(r"Paired lensing response: $C_L^{q\phi}$", loc="left")
    axes[0].set_ylabel(r"$10^7 L(L+1)\Delta C_L^{q\phi}$")
    axes[1].set_title("Paired marked phase-control residual", loc="left")
    axes[1].set_ylabel(r"$10^7\Delta C_L$")
    axes[2].set_title("Paired normalized connected-q component", loc="left")
    axes[2].set_ylabel(r"$10^7\Delta C_L$")
    axes[2].set_xlabel(r"Output multipole $L$")
    for ax in axes:
        ax.axhline(0, color="0.3", lw=0.8)
        ax.axvline(200, color="0.45", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlim(40, 3000)
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, ncol=2)
    fig.suptitle("Same CMB seed and output bins; two input temperature cuts", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "temperature_cut_comparison.png", dpi=190)
    plt.close(fig)
    with (OUT / "temperature_cut_comparison.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("ell_min", "ell_max", *(
            f"lmin_{cut}_{name}" for cut in FOLDERS for name in results[cut]
        )))
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            writer.writerow((lo, hi - 1, *(
                results[cut][name][i] for cut in FOLDERS for name in results[cut]
            )))
    print(f"Saved temperature-cut comparison in {OUT}")


if __name__ == "__main__":
    main()
