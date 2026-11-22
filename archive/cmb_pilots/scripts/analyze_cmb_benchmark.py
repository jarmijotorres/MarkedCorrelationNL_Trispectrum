#!/usr/bin/env python3
"""Summarize the paired p=2 CMB benchmark without treating surrogates as seeds."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


FIELDS = ("marked_minus_phase_controls", "connected_q_minus_phase_controls",
          "connected_q_component_minus_phase_controls")
KINDS = ("lensed", "unlensed")


def mean_sem(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return values.mean(axis=0), values.std(axis=0, ddof=1) / np.sqrt(values.shape[0])


def load_data(path: Path, nsims: int, nbins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return [pair, kind, bin, metric], ell bounds, and seed by pair."""
    metrics = (*FIELDS, "temperature_cl")
    values = np.full((nsims, 2, nbins, len(metrics)), np.nan)
    bounds = np.full((nbins, 2), -1, dtype=int)
    seeds = np.full(nsims, -1, dtype=int)
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            pair, band = int(row["pair"]), int(row["bin"])
            kind = KINDS.index(row["kind"])
            if not (0 <= pair < nsims and 0 <= band < nbins):
                raise ValueError("Pair or bin index lies outside benchmark metadata")
            if np.isfinite(values[pair, kind, band, 0]):
                raise ValueError(f"Duplicate pair/kind/bin: {pair}/{kind}/{band}")
            values[pair, kind, band] = [float(row[name]) for name in metrics]
            current = (int(row["ell_min"]), int(row["ell_max"]))
            if bounds[band, 0] >= 0 and tuple(bounds[band]) != current:
                raise ValueError(f"Inconsistent ell bounds in bin {band}")
            bounds[band] = current
            seed = int(row["seed"])
            if seeds[pair] >= 0 and seeds[pair] != seed:
                raise ValueError(f"Inconsistent seed for pair {pair}")
            seeds[pair] = seed
    if not np.isfinite(values).all() or (bounds < 0).any() or (seeds < 0).any():
        raise ValueError("Benchmark CSV has missing or non-finite measurements")
    if len(set(seeds)) != nsims:
        raise ValueError("Benchmark seeds must be independent across pairs")
    if (values[..., 3] <= 0).any():
        raise ValueError("Temperature bandpowers must be positive for normalization")
    return values, bounds, seeds


def summarize(values: np.ndarray, bounds: np.ndarray, output: Path) -> None:
    # Divide each seed by its own C_ell^TT. This is a plotting normalization,
    # while the raw C_ell estimates and paired differences remain in the CSV.
    temp = values[..., 3]
    pair_temp = temp.mean(axis=1)
    headings = ["bin", "ell_min", "ell_max", "descriptor", "lensed_mean", "lensed_sem",
                "unlensed_mean", "unlensed_sem", "paired_mean", "paired_sem",
                "lensed_over_tt_mean", "lensed_over_tt_sem",
                "unlensed_over_tt_mean", "unlensed_over_tt_sem",
                "paired_over_tt_mean", "paired_over_tt_sem"]
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headings)
        for band, (lo, hi) in enumerate(bounds):
            for name, data in (
                ("marked_residual", values[:, :, band, 0]),
                ("connected_q_residual", values[:, :, band, 1]),
                ("normalized_q_connected_component", values[:, :, band, 2]),
                ("marked_remainder", values[:, :, band, 0] - values[:, :, band, 2]),
            ):
                paired = data[:, 0] - data[:, 1]
                lens, lens_sem = mean_sem(data[:, 0])
                null, null_sem = mean_sem(data[:, 1])
                delta, delta_sem = mean_sem(paired)
                lens_n, lens_n_sem = mean_sem(data[:, 0] / temp[:, 0, band])
                null_n, null_n_sem = mean_sem(data[:, 1] / temp[:, 1, band])
                delta_n, delta_n_sem = mean_sem(paired / pair_temp[:, band])
                writer.writerow([band, lo, hi, name, lens, lens_sem, null, null_sem,
                                 delta, delta_sem, lens_n, lens_n_sem, null_n, null_n_sem,
                                 delta_n, delta_n_sem])


def plot(values: np.ndarray, bounds: np.ndarray, output: Path) -> None:
    centers = bounds.mean(axis=1)
    temp = values[..., 3]
    pair_temp = temp.mean(axis=1)
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 9.2), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 1.15]})
    colors = ("#bb4b2f", "#24718c")
    for ax, field_index, title in zip(
        axes[:2], (0, 1),
        (r"Marked: $\Delta C_\ell^{\delta_m\delta_m}$",
         r"Quadratic: $\Delta(C_\ell^{qq}-W_\ell)$")
    ):
        for kind_index, kind in enumerate(KINDS):
            normalized = values[:, kind_index, :, field_index] / temp[:, kind_index, :]
            mean = normalized.mean(axis=0)
            ax.plot(centers, mean, marker=("o", "s")[kind_index],
                    ms=4, lw=1.2, color=colors[kind_index], label=kind)
        ax.axhline(0, color="0.3", lw=0.8)
        ax.set_ylabel(r"Residual / $C_\ell^{\delta_T\delta_T}$")
        ax.set_title(title, loc="left", fontsize=11)
        ax.grid(alpha=0.18)
    axes[0].legend(frameon=False, ncol=2, loc="lower right")
    # Use the marked panel's scale for the connected projection comparison.
    axes[1].set_ylim(axes[0].get_ylim())
    difference = values[:, 0, :, :3] - values[:, 1, :, :3]
    for label, data, color, marker in (
        ("Marked", difference[..., 0], "#202c40", "o"),
        (r"Normalized $q$ connected part", difference[..., 2], "#9d659f", "s"),
        ("Remainder", difference[..., 0] - difference[..., 2], "#42a07a", "^"),
    ):
        mean = (data / pair_temp).mean(axis=0)
        axes[2].plot(centers, mean, marker=marker, ms=4, lw=1.2,
                     color=color, label=label)
    axes[2].axhline(0, color="0.3", lw=0.8)
    axes[2].set_ylabel(r"Paired difference / $C_\ell^{\delta_T\delta_T}$")
    axes[2].set_title("Lensed minus unlensed, after each phase control", loc="left", fontsize=11)
    axes[2].legend(frameon=False, ncol=3, fontsize=8, loc="upper right")
    axes[2].grid(alpha=0.18)
    axes[2].set_xlabel(r"Multipole $\ell$ (bin centers)")
    axes[2].set_xlim(bounds[0, 0] - 6, bounds[-1, 1] + 6)
    fig.suptitle(f"CMB p=2 mark: {values.shape[0]} paired seeds, Nside=128",
                 fontsize=13)
    fig.text(0.5, 0.01, f"Points are means over {values.shape[0]} paired seeds, plotted at the true bin centers.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 0.965))
    fig.savefig(output, dpi=190)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("archive/cmb_pilots/run_cmb_benchmark_128_p2"))
    args = parser.parse_args()
    metadata = json.loads((args.input_dir / "benchmark.json").read_text())
    values, bounds, _ = load_data(args.input_dir / "bandpowers.csv",
                                  metadata["nsims"], metadata["nbins"])
    if metadata["nsims"] < 2:
        parser.error("At least two independent pairs are needed for a seed-level standard error")
    if metadata.get("mark_exponent") != 2:
        parser.error("This analyzer expects the p=2 density-style CMB mark")
    summarize(values, bounds, args.input_dir / "summary_by_bin.csv")
    plot(values, bounds, args.input_dir / "benchmark_results.png")
    print(f"Analyzed {values.shape[0]} independent pairs in {values.shape[2]} ell bins")
    print(f"Wrote {args.input_dir / 'summary_by_bin.csv'}")
    print(f"Wrote {args.input_dir / 'benchmark_results.png'}")


if __name__ == "__main__":
    main()
