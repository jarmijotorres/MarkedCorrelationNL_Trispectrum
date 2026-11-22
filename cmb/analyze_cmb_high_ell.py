#!/usr/bin/env python3
"""Inspect the full output-L range of one paired ACT-like CMB realization."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cmb.benchmark_cmb_ensemble import bin_spectrum


SPECTRA = ("marked_residual", "connected_q_residual", "q_component_residual")
COLORS = {"lensed": "#bb4b2f", "unlensed": "#24718c"}


def correlation(cl: np.ndarray, theta_arcmin: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Full-sky isotropic Legendre transform, retaining lo <= L < hi."""
    ell = np.arange(lo, hi)
    coefficients = np.zeros(hi)
    coefficients[lo:hi] = (2 * ell + 1) * cl[lo:hi] / (4 * np.pi)
    return np.polynomial.legendre.legval(np.cos(np.deg2rad(theta_arcmin / 60)), coefficients)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, default=Path("run_cmb_angular_consistency_2048"))
    folder = parser.parse_args().folder
    meta = json.loads((folder / "benchmark.json").read_text())
    if meta["nsims"] != 1:
        raise ValueError("This diagnostic expects one paired realization")
    edges = np.asarray(meta["bin_edges"], dtype=int)
    centers = (edges[:-1] + edges[1:] - 1) / 2
    with np.load(folder / "unbinned_spectra.npz") as saved:
        spectra = {key: saved[key] for key in saved.files}
    data = {
        kind: {name: bin_spectrum(spectra[f"pair_0_{kind}_{name}"], edges)
               for name in (*SPECTRA, "q_phi")}
        for kind in COLORS
    }
    pair = {name: data["lensed"][name] - data["unlensed"][name] for name in SPECTRA}

    fig, axes = plt.subplots(4, 1, figsize=(9.3, 11.5), sharex=True)
    titles = (
        r"Marked phase-control residual $\Delta C_L^{\delta_m\delta_m}$",
        r"Quadratic phase-control residual $\Delta(C_L^{qq}-W_L)$",
    )
    for index, (name, title) in enumerate(zip(SPECTRA, titles)):
        ax = axes[index]
        for kind, color in COLORS.items():
            ax.plot(centers, data[kind][name] * 1e7, color=color, lw=1.2,
                    marker=".", ms=3, label=kind)
        ax.set_title(title, loc="left")
        ax.set_ylabel(r"$10^7 C_L$")
        ax.legend(frameon=False, ncol=2)
    axes[2].plot(centers, pair["marked_residual"] * 1e7,
                 color="#202c40", lw=1.25, label="marked")
    axes[2].plot(centers, pair["q_component_residual"] * 1e7,
                 color="#9d659f", lw=1.25, label="normalized connected q")
    axes[2].plot(centers, (pair["marked_residual"] - pair["q_component_residual"]) * 1e7,
                 color="#42a07a", lw=1.25, label="remainder")
    axes[2].set_ylabel(r"Paired $10^7 C_L$")
    axes[2].set_title("Lensed minus unlensed, after phase controls", loc="left")
    axes[2].legend(frameon=False, ncol=3, fontsize=8)
    for kind, color in COLORS.items():
        axes[3].plot(centers, 1e7 * centers * (centers + 1) * data[kind]["q_phi"],
                     color=color, lw=1.2, marker=".", ms=3, label=kind)
    axes[3].set_ylabel(r"$10^7 L(L+1)C_L^{q\phi}$")
    axes[3].set_title("Quadratic field versus known lensing potential", loc="left")
    axes[3].set_xlabel(r"Output multipole $L$")
    axes[3].legend(frameon=False, ncol=2)
    for ax in axes:
        ax.axhline(0, color="0.3", lw=0.8)
        ax.axvline(200, color="0.45", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.grid(alpha=0.15)
    axes[3].set_xlim(edges[0], edges[-1])
    fig.suptitle(
        fr"CMB mark: input temperature $\ell\geq{meta['temperature_lmin']}$, one pair",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(folder / "high_ell_results.png", dpi=190)
    plt.close(fig)

    theta = np.geomspace(2, 180, 400)
    ranges = {"low": (40, 200), "wide": (40, int(edges[-1]))}
    xi = {
        label: {
            kind: {name: correlation(spectra[f"pair_0_{kind}_{name}"], theta, *limits)
                   for name in SPECTRA}
            for kind in COLORS
        }
        for label, limits in ranges.items()
    }
    fig, axes = plt.subplots(3, 1, figsize=(9.3, 9.5), sharex=True)
    for ax, name, title in zip(
        axes,
        ("marked_residual", "connected_q_residual", "q_component_residual"),
        ("Marked residual", "Connected q residual", "Paired normalized q component"),
    ):
        if name == "q_component_residual":
            for label, linestyle in (("low", "--"), ("wide", "-")):
                values = xi[label]["lensed"][name] - xi[label]["unlensed"][name]
                ax.plot(theta, values * 1e4, color="#9d659f", ls=linestyle,
                        lw=1.4, label=f"{label} L range")
        else:
            for kind, color in COLORS.items():
                for label, linestyle in (("low", "--"), ("wide", "-")):
                    ax.plot(theta, xi[label][kind][name] * 1e4, color=color,
                            ls=linestyle, lw=1.4, label=f"{kind}, {label} L range")
        ax.axhline(0, color="0.3", lw=0.8)
        ax.set_title(title, loc="left")
        ax.set_ylabel(r"$10^4\,\xi(\theta)$")
        ax.set_xscale("log")
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8, ncol=2)
    axes[-1].set_xlabel(r"Angular separation $\theta$ [arcmin]")
    fig.suptitle("Band-limited correlations: dashed 40–199, solid 40–2999", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(folder / "angular_correlations_extended.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(9.3, 9), sharex=True)
    paired_curves = (
        ("marked_residual", "Paired marked residual", "#202c40"),
        ("connected_q_residual", "Paired connected q residual", "#bb4b2f"),
        ("q_component_residual", "Paired normalized q component", "#9d659f"),
    )
    for ax, (name, title, color) in zip(axes, paired_curves):
        for label, linestyle in (("low", "--"), ("wide", "-")):
            values = xi[label]["lensed"][name] - xi[label]["unlensed"][name]
            ax.plot(theta, values * 1e4, color=color, ls=linestyle, lw=1.5,
                    label=f"{label} L range")
        ax.axhline(0, color="0.3", lw=0.8)
        ax.set_title(title, loc="left")
        ax.set_ylabel(r"Paired $10^4\,\xi(\theta)$")
        ax.set_xscale("log")
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, ncol=2)
    axes[-1].set_xlabel(r"Angular separation $\theta$ [arcmin]")
    fig.suptitle("Lensed minus unlensed angular response after phase controls", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(folder / "paired_angular_response.png", dpi=190)
    plt.close(fig)

    with (folder / "angular_correlations_extended.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("theta_arcmin", *(
            f"{label}_{kind}_{name}"
            for label in ranges for kind in COLORS for name in SPECTRA
        )))
        for i, angle in enumerate(theta):
            writer.writerow((angle, *(
                xi[label][kind][name][i]
                for label in ranges for kind in COLORS for name in SPECTRA
            )))
    print(f"Saved extended-L figures and correlation table in {folder}")


if __name__ == "__main__":
    main()
