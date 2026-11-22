#!/usr/bin/env python3
"""Rebuild the CMB connected-q residual by subtracting Wick terms in angle space."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import roots_legendre

from cmb.analyze_cmb_high_ell import correlation
from cmb.benchmark_cmb_ensemble import bin_spectrum
from cmb.core import correlation_to_cl


KINDS = ("lensed", "unlensed")
PROJECTIONS = ("connected_q_residual", "q_component_residual")


def angular_terms(saved, prefix, mu, window):
    """Return the angular q correlation and disconnected four-point terms."""
    tt = saved[f"{prefix}_tt"]
    qq = saved[f"{prefix}_qq"]
    factors = (2 * np.arange(tt.size) + 1) / (4 * np.pi)
    tt_coeff = factors * tt
    xi_tt = np.polynomial.legendre.legval(mu, tt_coeff)
    xi_ts = np.polynomial.legendre.legval(mu, tt_coeff * window)
    xi_ss = np.polynomial.legendre.legval(mu, tt_coeff * window**2)
    xi_qq = np.polynomial.legendre.legval(mu, factors * qq)
    return xi_qq, xi_tt * xi_ss + xi_ts**2, float(saved[f"{prefix}_q_coefficient"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    folder = args.folder
    meta = json.loads((folder / "benchmark.json").read_text())
    if meta["nsims"] != 1:
        parser.error("the comparison figure expects a single paired realization")
    edges = np.asarray(meta["bin_edges"], dtype=int)
    lo, hi = int(edges[0]), int(edges[-1])
    if not lo < 200 < hi:
        parser.error("expected a measured range spanning output L=200")
    ell = np.arange(meta["lmax"] + 1)
    window = np.exp(-0.5 * ell * (ell + 1) * np.deg2rad(meta["smoothing_degrees"])**2)
    # The Wick product can have degree 2*lmax. Integrate its product with
    # every retained P_L over the entire mu interval, including unplotted angles.
    mu, weights = roots_legendre(meta["lmax"] + (hi + 1) // 2)
    saved_cl = {}
    angle_cl = {}
    stats = {}
    with np.load(folder / "unbinned_spectra.npz") as saved:
        for kind in KINDS:
            prefix = f"pair_0_{kind}"
            qq, wick, kq = angular_terms(saved, prefix, mu, window)
            controls = [angular_terms(saved, f"{prefix}_control_{i}", mu, window)
                        for i in range(meta["surrogates"])]
            raw = qq - np.mean([cq for cq, _, _ in controls], axis=0)
            wick_residual = wick - np.mean([cw for _, cw, _ in controls], axis=0)
            connected = raw - wick_residual
            normalized = kq * (qq - wick) - np.mean(
                [ck * (cq - cw) for cq, cw, ck in controls], axis=0)
            for name, curve in zip(PROJECTIONS, (connected, normalized)):
                saved_cl[kind, name] = saved[f"{prefix}_{name}"]
                angle_cl[kind, name] = correlation_to_cl(curve, mu, weights, hi - 1)
                band = slice(lo, hi)
                error = np.max(np.abs(angle_cl[kind, name][band] - saved_cl[kind, name][band]))
                peak = np.max(np.abs(saved_cl[kind, name][band]))
                stats[f"{kind}_{name}_max_abs_cl_difference"] = float(error)
                np.testing.assert_allclose(angle_cl[kind, name][band], saved_cl[kind, name][band],
                                           atol=max(peak * 1e-6, 1e-15), rtol=0)
            stats[f"{kind}_max_abs_angle_wick_residual"] = float(np.max(np.abs(wick_residual)))

    theta = np.geomspace(2, 180, 450)
    q, kq = PROJECTIONS
    harmonic = saved_cl["lensed", q] - saved_cl["unlensed", q]
    angular = angle_cl["lensed", q] - angle_cl["unlensed", q]
    full = correlation(harmonic, theta, lo, hi)
    rebuilt = correlation(angular, theta, lo, hi)
    low = correlation(harmonic, theta, lo, 200)
    high = correlation(harmonic, theta, 200, hi)
    stats["paired_max_abs_xi_difference"] = float(np.max(np.abs(full - rebuilt)))
    for angle in (2, 5, 10, 15, 20, 30, 60):
        sample = np.array([angle], dtype=float)
        stats[f"paired_xi_at_{angle}_arcmin"] = {
            "low": float(correlation(harmonic, sample, lo, 200)[0]),
            "high": float(correlation(harmonic, sample, 200, hi)[0]),
            "total": float(correlation(harmonic, sample, lo, hi)[0]),
        }

    fig, axes = plt.subplots(3, 1, figsize=(8.7, 9.4))
    centers = (edges[:-1] + edges[1:] - 1) / 2
    axes[0].semilogx(centers, 1e7 * bin_spectrum(harmonic, edges), color="#202c40",
                     lw=1.5, label="Harmonic calculation")
    axes[0].semilogx(centers, 1e7 * bin_spectrum(angular, edges), "+", color="#bb4b2f",
                     ms=4, label="Recovered by integrating the angular correlation")
    axes[0].set(xlim=(lo, hi), xlabel=r"Output multipole $L$", ylabel=r"$10^7\,\Delta T_{q,L}$",
                title="Paired four-point residual in harmonic space")
    axes[1].semilogx(theta, 1e4 * low, "--", color="#24718c", label=f"$L={lo}$–199")
    axes[1].semilogx(theta, 1e4 * high, ":", color="#9d659f", lw=2, label=f"$L=200$–{hi-1}")
    axes[1].semilogx(theta, 1e4 * full, color="#202c40", label="Sum")
    axes[1].semilogx(theta[::24], 1e4 * rebuilt[::24], "+", color="#bb4b2f", ms=5,
                     label="Angle calculation with the same output band")
    axes[1].set(ylabel=r"$10^4\,\Delta T_q(\theta)$", title="How output bands add in angle space")

    for kind, color in zip(KINDS, ("#bb4b2f", "#24718c")):
        curve = correlation(saved_cl[kind, kq], theta, lo, hi)
        axes[2].semilogx(theta, 1e4 * curve, color=color, label=kind)
    marked_q = saved_cl["lensed", kq] - saved_cl["unlensed", kq]
    axes[2].semilogx(theta, 1e4 * correlation(marked_q, theta, lo, hi), color="#202c40",
                     ls="--", label="Lensed minus unlensed")
    axes[2].set(xlabel=r"Angular separation $\theta$ [arcmin]", ylabel=r"$10^4\,\Delta(K_q T_q)$",
                title="Connected q contribution to the marked correlation")
    for ax in axes:
        ax.axhline(0, color="0.5", lw=0.8)
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=9)
    for ax in axes[1:]:
        ax.axvline(30, color="0.55", lw=0.8, alpha=0.7)
        ax.set_xlim(theta[0], theta[-1])
    axes[1].set_xlabel(r"Angular separation $\theta$ [arcmin]")
    fig.suptitle(r"Same lensed/unlensed pair and phase controls in both bases", fontsize=13)
    fig.tight_layout()
    fig.savefig(folder / "angular_projection_check.png", dpi=190)
    plt.close(fig)
    (folder / "angular_projection_check.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
