#!/usr/bin/env python3
"""Analyse saved marked-spectrum ensembles without generating new maps."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from cmb.core import bin_spectrum


def analyze(folder: Path) -> str:
    """Write the decomposition figure and report for a completed ensemble."""
    import matplotlib.pyplot as plt

    folder = Path(folder)
    config = json.loads((folder / "config_snapshot.json").read_text())
    edges = np.asarray(config["bin_edges"])
    centers = bin_spectrum(np.arange(config["lmax"] + 1), edges)
    full, quadratic, unlensed = [], [], []
    for pair in sorted(folder.glob("pair_*/unbinned_spectra.npz")):
        with np.load(pair) as saved:
            full.append(bin_spectrum(
                saved["pair_0_lensed_marked_residual"]
                - saved["pair_0_unlensed_marked_residual"], edges))
            quadratic.append(bin_spectrum(
                saved["pair_0_lensed_q_component_residual"]
                - saved["pair_0_unlensed_q_component_residual"], edges))
            unlensed.append(bin_spectrum(saved["pair_0_unlensed_marked_residual"],
                                          edges))
    full, quadratic, unlensed = map(np.asarray, (full, quadratic, unlensed))
    if len(full) < 2:
        raise ValueError("Need at least two paired skies")
    remainder = full - quadratic

    fig, axes = plt.subplots(2, 1, figsize=(9, 8),
                             constrained_layout=True)
    colors = (("#202c40", full, "full marked"),
              ("#9d659f", quadratic, "connected q contribution"),
              ("#42a07a", remainder, "remainder"))
    for color, values, label in colors:
        axes[0].plot(centers, values.mean(axis=0) * 1e7, lw=2,
                     color=color, label=label)
    axes[0].set(xlabel="Output multipole L", ylabel=r"$10^7\Delta C_L$",
                title="Paired marked response and its connected-q component")
    for color, values, label in colors:
        low, high = (100, 600)
        chosen = (edges[:-1] >= low) & (edges[1:] <= high)
        # Plot paired values in the primary band for each independent sky.
        weights = edges[1:][chosen] ** 2 - edges[:-1][chosen] ** 2
        band = np.average(values[:, chosen], axis=1, weights=weights)
        axes[1].plot(np.arange(1, len(values) + 1), band * 1e7,
                     "o-", color=color, label=label)
    axes[1].set(xlabel="Paired sky", ylabel=r"$10^7$ bandpower",
                title=r"Independent skies, $100\leq L<600$")
    axes[1].set_xticks(np.arange(1, len(full) + 1))
    for ax in axes:
        ax.axhline(0, color="0.45", lw=0.8)
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
    fig.savefig(folder / "marked_decomposition.png", dpi=180)
    plt.close(fig)

    bands = [(100, 600), (600, 1200), (640, 800), (800, 1000),
             (1000, 1200)]
    report = [f"# Saved-sky diagnostic ({len(full)} pairs)", ""]
    report.append("Values are phase-control-subtracted lensed minus unlensed "
                  "unless identified otherwise. Sample SD measures sky-to-sky "
                  "spread; SEM is SD divided by the square root of the number "
                  "of paired skies.")
    report.extend(["", "| Output L | Marked mean | Sample SD | SEM | "
                   "Connected q mean | Remainder mean | Positive pairs |",
                   "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
    paired_bands = []
    for low, high in bands:
        selected = (edges[:-1] >= low) & (edges[1:] <= high)
        weights = edges[1:][selected] ** 2 - edges[:-1][selected] ** 2
        mark = np.average(full[:, selected], axis=1, weights=weights)
        quad = np.average(quadratic[:, selected], axis=1, weights=weights)
        rem = mark - quad
        paired_bands.append(mark)
        report.append(f"| {low}–{high} | {mark.mean():.3e} | "
                      f"{mark.std(ddof=1):.3e} | "
                      f"{mark.std(ddof=1) / np.sqrt(len(mark)):.3e} | "
                      f"{quad.mean():.3e} | {rem.mean():.3e} | "
                      f"{np.count_nonzero(mark > 0)}/{len(mark)} |")
    broad = np.array(paired_bands[:2]).T
    covariance = np.cov(broad, rowvar=False)
    correlation = np.corrcoef(broad, rowvar=False)
    report.extend(["", "## Broad-band sample covariance", "",
                   "Rows and columns are 100–600, then 600–1200. "
                   "Covariance units are C_L squared.",
                   "",
                   f"- Covariance: {covariance.tolist()}",
                   f"- Correlation: {correlation[0, 1]:.3f}",
                   f"- With {len(full)} skies and {len(centers)} bins, the "
                   "many-bin covariance has rank at most the smaller of "
                   f"{len(full) - 1} and {len(centers)}."])
    vector_file = folder / "ensemble_vectors.npz"
    if vector_file.exists():
        with np.load(vector_file) as vectors:
            if "paired_null_covariance" in vectors:
                null_covariance = vectors["paired_null_covariance"]
                report.extend([
                    "", "## Paired phase-randomized null", "",
                    "For each sky, matched lensed-minus-unlensed phase-control "
                    "contrasts estimate the Gaussian-phase covariance. The "
                    "variance of one target contrast minus the mean of K "
                    "controls is estimated as (1+1/K) times that covariance. "
                    "The null-mean SEM divides the single-sky SD by the square "
                    "root of the number of independent skies. This is a "
                    "conditional phase null, not an unlensed/no-lensing "
                    "simulation ensemble.",
                    "",
                    "| Output L | Single-sky null SD | Null-mean SEM | "
                    "Observed paired mean |",
                    "| --- | ---: | ---: | ---: |",
                ])
                for (low, high), observed in zip(bands[:2], paired_bands[:2]):
                    selected = (edges[:-1] >= low) & (edges[1:] <= high)
                    weights = (edges[1:][selected] ** 2
                               - edges[:-1][selected] ** 2).astype(float)
                    weights /= weights.sum()
                    covariance = null_covariance[np.ix_(selected, selected)]
                    null_sd = float(np.sqrt(weights @ covariance @ weights))
                    report.append(
                        f"| {low}–{high} | {null_sd:.3e} | "
                        f"{null_sd / np.sqrt(len(full)):.3e} | "
                        f"{observed.mean():.3e} |"
                    )
    bump = paired_bands[3] - 0.5 * (paired_bands[2] + paired_bands[4])
    report.extend(["", "## High-L shape check", "",
                   "Define the 800–1000 band minus the average of 640–800 "
                   "and 1000–1200. This is a descriptive shape check, "
                   "chosen from the ten-sky pilot; independent future skies "
                   "can test it.",
                   f"- Mean contrast: {bump.mean():.3e}; "
                   f"sample SD: {bump.std(ddof=1):.3e}; "
                   f"positive pairs: {np.count_nonzero(bump > 0)}/{len(bump)}."])
    unlensed_band = []
    for low, high in bands[:2]:
        selected = (edges[:-1] >= low) & (edges[1:] <= high)
        unlensed_band.append(np.average(unlensed[:, selected], axis=1,
                                        weights=(edges[1:][selected] ** 2
                                                 - edges[:-1][selected] ** 2)))
    report.extend(["", "## Unlensed phase-control residual", ""])
    for (low, high), values in zip(bands[:2], unlensed_band):
        report.append(f"- {low}–{high}: mean {values.mean():.3e}, "
                      f"sample SD {values.std(ddof=1):.3e}.")
    result = "\n".join(report) + "\n"
    (folder / "saved_sky_analysis.md").write_text(result)
    return result


def main() -> None:
    import matplotlib
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", nargs="?", type=Path,
                        default=Path("run_cmb_marked_10_2048"))
    args = parser.parse_args()
    print(analyze(args.folder), end="")


if __name__ == "__main__":
    main()
