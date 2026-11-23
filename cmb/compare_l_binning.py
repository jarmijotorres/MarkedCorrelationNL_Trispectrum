"""Compare original and grouped L bins in the unlensed single-sky null test.

Usage: python -m cmb.compare_l_binning PATH/ensemble_vectors.npz
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


GROUPINGS = {
    "Original 43 bins": None,
    "Log-like 8 bins": np.array([40, 80, 120, 180, 260, 380, 560, 800, 1200]),
    "High L: 15 bins": np.arange(600, 1201, 40),
    "High L: 3 bins": np.array([600, 800, 1000, 1200]),
}


def group_spectra(spectra, original_edges, new_edges):
    if new_edges is None:
        return spectra
    grouped = []
    for lo, hi in zip(new_edges[:-1], new_edges[1:]):
        indices = np.flatnonzero((original_edges[:-1] >= lo) & (original_edges[1:] <= hi))
        if not len(indices) or original_edges[indices[0]] != lo or original_edges[indices[-1] + 1] != hi:
            raise ValueError(f"[{lo}, {hi}) does not align with saved bin edges")
        # Full-sky C_L bandpowers use sum_{L=lo}^{hi-1}(2L+1) = hi^2-lo^2.
        weights = (original_edges[indices + 1] ** 2 - original_edges[indices] ** 2).astype(float)
        grouped.append(spectra[:, indices] @ (weights / weights.sum()))
    return np.stack(grouped, axis=1)


def crossfit_scores(unlensed, lensed, folds):
    n = len(unlensed)
    null_scores, lens_scores = np.empty(n), np.empty(n)
    for test in folds:
        train = np.setdiff1d(np.arange(n), test)
        mean = unlensed[train].mean(axis=0)
        covariance = np.atleast_2d(np.cov(unlensed[train], rowvar=False))
        covariance = 0.9 * covariance + 0.1 * np.diag(np.diag(covariance))
        for data, scores in ((unlensed, null_scores), (lensed, lens_scores)):
            residual = data[test] - mean
            scores[test] = np.einsum(
                "ni,ni->n", residual, np.linalg.solve(covariance, residual.T).T
            )
    return null_scores, lens_scores


def run(vector_path: Path, output_path: Path):
    with np.load(vector_path) as data:
        edges = data["bin_edges"]
        unlensed = data["unlensed"]
        lensed = data["lensed"]
    n = len(unlensed)
    folds = np.array_split(np.random.default_rng(42).permutation(n), 5)
    results = {}
    for name, grouping in GROUPINGS.items():
        null, lens = crossfit_scores(group_spectra(unlensed, edges, grouping),
                                     group_spectra(lensed, edges, grouping), folds)
        threshold = np.percentile(null, 95)
        results[name] = (null, lens, threshold)
        print(f"{name}: {(lens > threshold).sum()}/{n} lensed skies above null 95th percentile")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, (name, (null, lens, threshold)) in zip(axes.flat, results.items()):
        bins = np.histogram_bin_edges(np.r_[null, lens], bins=18)
        ax.hist(null, bins=bins, alpha=0.55, color="tab:blue", label="Unlensed")
        ax.hist(lens, bins=bins, alpha=0.55, color="tab:orange", label="Lensed")
        ax.axvline(threshold, color="tab:blue", ls="--", lw=1.5)
        ax.set_title(f"{name}: {(lens > threshold).sum()}/{n} above null 95th percentile")
        ax.set_xlabel("Covariance-weighted distance from unlensed mean")
        ax.set_ylabel("Number of skies")
        ax.legend(fontsize=8)
    fig.suptitle("Exploratory L-binning comparison · same 5 folds · 10% diagonal covariance shrinkage")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    np.savez_compressed(output_path.with_suffix(".npz"),
                        log8_edges=GROUPINGS["Log-like 8 bins"],
                        **{f"{i}_{kind}": values[k]
                           for i, values in enumerate(results.values())
                           for kind, k in (("null", 0), ("lensed", 1), ("threshold", 2))})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vectors", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run(args.vectors, args.output or args.vectors.with_name("l_binning_null_comparison.png"))
