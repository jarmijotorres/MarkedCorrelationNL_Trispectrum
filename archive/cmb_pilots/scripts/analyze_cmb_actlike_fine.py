#!/usr/bin/env python3
"""Plot fine multipole bins and their band-limited angular correlations."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cmb.benchmark_cmb_ensemble import bin_spectrum


FOLDER = Path("archive/cmb_pilots/run_cmb_actlike_fine_2048")


def legendre_correlation(cl: np.ndarray, theta_degrees: np.ndarray,
                         ell_min: int, ell_max: int) -> np.ndarray:
    """Transform only the retained multipoles to an isotropic xi(theta)."""
    coeff = np.zeros(ell_max + 1)
    ell = np.arange(ell_min, ell_max + 1)
    coeff[ell] = (2 * ell + 1) * cl[ell] / (4 * np.pi)
    return np.polynomial.legendre.legval(np.cos(np.deg2rad(theta_degrees)), coeff)


def main() -> None:
    meta = json.loads((FOLDER / "benchmark.json").read_text())
    if meta["nsims"] != 1:
        raise ValueError("expected one paired realization")
    edges = np.asarray(meta["bin_edges"])
    if not np.all(np.diff(edges) == 10):
        raise ValueError("expected bins of width Delta ell = 10")
    spectra = np.load(FOLDER / "unbinned_spectra.npz")
    names = ("marked_residual", "connected_q_residual", "q_component_residual")
    data = {
        kind: np.column_stack([
            bin_spectrum(spectra[f"pair_0_{kind}_{name}"], edges) for name in names
        ])
        for kind in ("lensed", "unlensed")
    }
    delta = data["lensed"] - data["unlensed"]
    centers = (edges[:-1] + edges[1:] - 1) / 2
    colors = ("#bb4b2f", "#24718c")
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for index, title in enumerate((r"Marked: $\Delta C_\ell^{\delta_m\delta_m}$",
                                   r"Quadratic: $\Delta(C_\ell^{qq}-W_\ell)$")):
        for kind, marker, color in zip(("lensed", "unlensed"), ("o", "s"), colors):
            axes[index].plot(centers, data[kind][:, index] * 1e7,
                             marker=marker, color=color, lw=1.2, ms=3.5,
                             label=kind)
        axes[index].axhline(0, color="0.3", lw=0.8)
        axes[index].set_title(title, loc="left", fontsize=11)
        axes[index].set_ylabel(r"Residual $C_\ell\times10^7$")
        axes[index].grid(alpha=0.17)
    limit = 1.1 * max(abs(data[kind][:, index] * 1e7).max()
                      for kind in data for index in (0, 1))
    axes[0].set_ylim(-limit, limit)
    axes[1].set_ylim(-limit, limit)
    axes[0].legend(frameon=False, ncol=2)
    for label, values, marker, color in (
        ("Marked", delta[:, 0], "o", "#202c40"),
        (r"Normalized $q$ connected part", delta[:, 2], "s", "#9d659f"),
        ("Remainder", delta[:, 0] - delta[:, 2], "^", "#42a07a"),
    ):
        axes[2].plot(centers, values * 1e7, marker=marker, color=color,
                     lw=1.2, ms=3.5, label=label)
    axes[2].axhline(0, color="0.3", lw=0.8)
    axes[2].set_title("Lensed minus unlensed, after each phase control", loc="left", fontsize=11)
    axes[2].set_ylabel(r"Paired difference $\times10^7$")
    axes[2].set_xlabel(r"Output multipole $\ell$ (bin width $\Delta\ell=10$)")
    axes[2].legend(frameon=False, fontsize=8, ncol=3)
    axes[2].grid(alpha=0.17)
    fig.suptitle("ACT-like CMB mark: one paired realization", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FOLDER / "fine_ell_results.png", dpi=190)
    plt.close(fig)

    theta = np.linspace(0.25, 10, 220)
    xi = {
        kind: np.column_stack([
            legendre_correlation(spectra[f"pair_0_{kind}_{name}"], theta,
                                 int(edges[0]), int(edges[-1]) - 1)
            for name in names
        ])
        for kind in ("lensed", "unlensed")
    }
    xi_delta = xi["lensed"] - xi["unlensed"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for index, title in enumerate((r"Marked: $\Delta\xi_{\delta_m\delta_m}(\theta)$",
                                   r"Quadratic: $\Delta T_q(\theta)$")):
        for kind, color in zip(("lensed", "unlensed"), colors):
            axes[index].plot(theta, xi[kind][:, index] * 1e4,
                             color=color, lw=1.5, label=kind)
        axes[index].axhline(0, color="0.3", lw=0.8)
        axes[index].set_title(title, loc="left", fontsize=11)
        axes[index].set_ylabel(r"Correlation $\times10^4$")
        axes[index].grid(alpha=0.17)
    limit = 1.1 * max(abs(xi[kind][:, index] * 1e4).max()
                      for kind in xi for index in (0, 1))
    axes[0].set_ylim(-limit, limit)
    axes[1].set_ylim(-limit, limit)
    axes[0].legend(frameon=False, ncol=2)
    for label, values, color in (
        ("Marked", xi_delta[:, 0], "#202c40"),
        (r"Normalized $q$ connected part", xi_delta[:, 2], "#9d659f"),
        ("Remainder", xi_delta[:, 0] - xi_delta[:, 2], "#42a07a"),
    ):
        axes[2].plot(theta, values * 1e4, color=color, lw=1.5, label=label)
    axes[2].axhline(0, color="0.3", lw=0.8)
    axes[2].set_title("Lensed minus unlensed, after each phase control", loc="left", fontsize=11)
    axes[2].set_ylabel(r"Paired difference $\times10^4$")
    axes[2].set_xlabel(r"Angular separation $\theta$ [degrees]")
    axes[2].legend(frameon=False, fontsize=8, ncol=3)
    axes[2].grid(alpha=0.17)
    fig.suptitle(r"Angular correlations from output multipoles $40\leq\ell<200$", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FOLDER / "angular_correlations.png", dpi=190)
    plt.close(fig)

    with (FOLDER / "angular_correlations.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("theta_degrees", "lensed_marked_residual", "unlensed_marked_residual",
                         "lensed_connected_q_residual", "unlensed_connected_q_residual",
                         "paired_marked", "paired_connected_q_component", "paired_remainder"))
        for i, angle in enumerate(theta):
            writer.writerow((angle, xi["lensed"][i, 0], xi["unlensed"][i, 0],
                             xi["lensed"][i, 1], xi["unlensed"][i, 1],
                             xi_delta[i, 0], xi_delta[i, 2], xi_delta[i, 0] - xi_delta[i, 2]))
    with (FOLDER / "bandpowers.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    qphi = {
        kind: np.array([float(row["q_phi_cl"]) for row in sorted(
            (row for row in rows if row["kind"] == kind), key=lambda row: int(row["bin"])
        )])
        for kind in ("lensed", "unlensed")
    }
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for kind, marker, color in zip(("lensed", "unlensed"), ("o", "s"), colors):
        ax.plot(centers, 1e7 * centers * (centers + 1) * qphi[kind],
                marker=marker, ms=4, lw=1.3, color=color, label=kind)
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xlabel(r"Output multipole $L$ (bin width $\Delta L=10$)")
    ax.set_ylabel(r"$10^7 L(L+1)C_L^{q\phi}$")
    ax.set_title("Quadratic field versus known input lensing potential", loc="left")
    ax.legend(frameon=False)
    ax.grid(alpha=0.17)
    fig.tight_layout()
    fig.savefig(FOLDER / "q_phi_fine.png", dpi=190)
    plt.close(fig)
    print(f"Wrote fine ell and angular correlation figures in {FOLDER}")


if __name__ == "__main__":
    main()
