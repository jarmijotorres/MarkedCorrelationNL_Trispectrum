"""Plot the physical unlensed null and the matched-pair response separately.

Usage: python -m cmb.plot_null_tests PATH/ensemble_vectors.npz
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_plot(vector_path: Path, output_path: Path) -> None:
    with np.load(vector_path) as data:
        centers = data["centers"]
        lensed = data["lensed"]
        unlensed = data["unlensed"]
        paired = data["paired_marked"]

    if lensed.shape != unlensed.shape or paired.shape != lensed.shape:
        raise ValueError("Lensed, unlensed, and paired vectors must have matching shapes")
    if not np.allclose(paired, lensed - unlensed, rtol=1e-9, atol=1e-15):
        raise ValueError("Stored paired vectors do not equal lensed minus unlensed")

    n = lensed.shape[0]
    null_mean = unlensed.mean(axis=0)
    null_lo, null_hi = np.percentile(unlensed, [2.5, 97.5], axis=0)
    lensed_mean = lensed.mean(axis=0)

    pair_mean = paired.mean(axis=0)
    pair_sem = paired.std(axis=0, ddof=1) / np.sqrt(n)

    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    scale = 1e8

    ax = axes[0]
    ax.fill_between(centers, scale * null_lo, scale * null_hi, color="tab:blue", alpha=0.22,
                    label="Unlensed null: central 95% of individual skies")
    ax.plot(centers, scale * null_mean, color="tab:blue", lw=1.7, label="Unlensed mean")
    ax.plot(centers, scale * lensed_mean, color="tab:orange", lw=2, label="Lensed mean")
    ax.set_title("1. Observable: lensed residual versus physical unlensed null")
    ax.set_ylabel(r"$\Delta C_L^{\rm mark}\;[10^{-8}]$")
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]
    ax.axhline(0, color="0.3", lw=1, label="Zero response")
    ax.fill_between(centers, scale * (pair_mean - pair_sem),
                    scale * (pair_mean + pair_sem), color="black", alpha=0.3,
                    label="Paired mean ± SEM")
    ax.plot(centers, scale * pair_mean, color="black", lw=2, label="Paired mean")
    ax.set_title("2. Simulation response: lensed residual minus matched unlensed residual")
    ax.set_ylabel(r"$D_L\;[10^{-8}]$")
    ax.set_xlabel(r"Marked-field multipole $L$")
    ax.legend(fontsize=8, loc="upper right")

    for ax in axes:
        ax.grid(alpha=0.2)
    axes[-1].set_xlim(centers.min(), centers.max())
    fig.suptitle(f"Marked-spectrum diagnostics from {n} matched CMB pairs", fontsize=14)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vectors", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    make_plot(args.vectors, args.output or args.vectors.with_name("marked_null_tests.png"))
