"""Exploratory single-sky null test using the unlensed marked-spectrum covariance.

Usage: python -m cmb.test_unlensed_null PATH/ensemble_vectors.npz
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def run(vector_path: Path, output_path: Path) -> None:
    with np.load(vector_path) as data:
        centers = data["centers"]
        edges = data["bin_edges"]
        unlensed = data["unlensed"]
        lensed = data["lensed"]
    if unlensed.shape != lensed.shape or unlensed.shape[1] != centers.size:
        raise ValueError("Incompatible lensed, unlensed, or bin arrays")

    n, d = unlensed.shape
    if n < d + 2:
        raise ValueError("More null maps than bins are needed for this covariance test")

    null_covariance = np.cov(unlensed, rowvar=False, ddof=1)
    sigma = np.sqrt(np.diag(null_covariance))
    correlation = null_covariance / np.outer(sigma, sigma)

    # Each fold uses only its other 80 unlensed skies to define the null.
    # Both held-out groups are scored with exactly the same fitted mean/covariance.
    rng = np.random.default_rng(42)
    folds = np.array_split(rng.permutation(n), 5)
    null_scores = np.empty(n)
    lensed_scores = np.empty(n)
    shrinkage = 0.10  # Stabilize inversion with 80 training maps and 43 bins.
    for test_idx in folds:
        train_idx = np.setdiff1d(np.arange(n), test_idx)
        train = unlensed[train_idx]
        mean = train.mean(axis=0)
        covariance = np.cov(train, rowvar=False, ddof=1)
        covariance = ((1 - shrinkage) * covariance
                      + shrinkage * np.diag(np.diag(covariance)))
        for maps, scores in ((unlensed, null_scores), (lensed, lensed_scores)):
            residual = maps[test_idx] - mean
            scores[test_idx] = np.einsum(
                "ni,ni->n", residual, np.linalg.solve(covariance, residual.T).T
            )

    # Ranking against held-out null scores is descriptive; the 100 skies and
    # shared fold fits do not support precise tail probabilities.
    empirical_p = (1 + (null_scores[None, :] >= lensed_scores[:, None]).sum(axis=1)) / (n + 1)
    null_self_p = (1 + (null_scores[None, :] >= null_scores[:, None]).sum(axis=1)) / (n + 1)
    threshold = np.percentile(null_scores, 95)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    im = ax0.pcolormesh(edges, edges, correlation, cmap="coolwarm", vmin=-1, vmax=1,
                        shading="flat")
    ax0.set_title(f"Unlensed null: {d}-bin correlation matrix")
    ax0.set_xlabel("Marked-field multipole $L$")
    ax0.set_ylabel("Marked-field multipole $L$")
    fig.colorbar(im, ax=ax0, label="Correlation")

    bins = np.histogram_bin_edges(np.r_[null_scores, lensed_scores], bins=20)
    ax1.hist(null_scores, bins=bins, alpha=0.55, label="Held-out unlensed skies", color="tab:blue")
    ax1.hist(lensed_scores, bins=bins, alpha=0.55, label="Held-out lensed skies", color="tab:orange")
    ax1.axvline(threshold, color="tab:blue", ls="--", lw=1.6,
                label="95th percentile of null scores")
    ax1.set_title("Single-sky distance from unlensed null")
    ax1.set_xlabel(r"$(X-\bar X_0)^T\widehat\Sigma_0^{-1}(X-\bar X_0)$")
    ax1.set_ylabel("Number of skies")
    ax1.legend(fontsize=8)
    fig.suptitle(f"Exploratory null test: {n} matched sky pairs; 5-fold validation; 10% diagonal covariance shrinkage")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    np.savez_compressed(output_path.with_suffix(".npz"), centers=centers,
                        null_covariance=null_covariance, null_correlation=correlation,
                        null_scores=null_scores, lensed_scores=lensed_scores,
                        lensed_empirical_p=empirical_p, null_empirical_p=null_self_p,
                        null_95th_score=threshold)
    print(f"Unlensed-null 95th-percentile score: {threshold:.2f}")
    print(f"Lensed skies above that threshold: {(lensed_scores > threshold).sum()}/{n}")
    print(f"Median lensed empirical null p: {np.median(empirical_p):.3f}")
    print(f"Saved {output_path} and {output_path.with_suffix('.npz')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vectors", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run(args.vectors, args.output or args.vectors.with_name("unlensed_covariance_null_test.png"))
