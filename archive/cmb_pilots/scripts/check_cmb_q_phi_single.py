#!/usr/bin/env python3
"""Cross-check an ACT-like quadratic temperature map against input lensing phi."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from lenspyx.lensing import synfast

from cmb.benchmark_cmb_ensemble import bin_spectrum


def main() -> None:
    folder = Path("archive/cmb_pilots/run_cmb_actlike_single_2048")
    config = json.loads((folder / "benchmark.json").read_text())
    if config["nsims"] != 1:
        raise ValueError("expected one seed pair")
    with (folder / "bandpowers.csv").open(newline="") as stream:
        seed = int(next(csv.DictReader(stream))["seed"])
    spectra = np.loadtxt(folder / "input_spectra.csv", delimiter=",", skiprows=1)
    tt, pp = spectra[:, 1], spectra[:, 2]
    nside, lmax = config["nside"], config["lmax"]
    ell = np.arange(lmax + 1)
    window = np.exp(-0.5 * ell * (ell + 1) * np.deg2rad(config["smoothing_degrees"]) ** 2)
    keep = (ell >= config["temperature_lmin"]).astype(float)
    maps, (alms, labels) = synfast(
        {"tt": tt, "pp": pp}, lmax=lmax,
        geometry=("healpix", {"nside": nside}), seed=seed,
        alm=True, nthreads=config["threads"], verbose=False,
    )
    phi_alm = alms[labels.index("p")]
    phi_cl = hp.alm2cl(phi_alm)
    unlensed = hp.alm2map(alms[labels.index("t")], nside, lmax=lmax)
    edges = np.asarray(config["bin_edges"])
    results = {}
    for kind, values in (("lensed", maps["T"]), ("unlensed", unlensed)):
        filtered_alm = hp.almxfl(hp.map2alm(values, lmax=lmax, iter=1), keep)
        t = hp.alm2map(filtered_alm, nside, lmax=lmax)
        del filtered_alm
        t -= t.mean()
        t /= t.std()
        t_alm = hp.map2alm(t, lmax=lmax, iter=3)
        s = hp.alm2map(hp.almxfl(t_alm, window), nside, lmax=lmax)
        q = t * s
        q -= q.mean()
        q_alm = hp.map2alm(q, lmax=lmax, iter=3)
        q_phi = hp.alm2cl(q_alm, phi_alm)
        q_q = hp.alm2cl(q_alm)
        results[kind] = (bin_spectrum(q_phi, edges), bin_spectrum(q_q, edges))
        del t, s, q, t_alm, q_alm
    phi_binned = bin_spectrum(phi_cl, edges)
    output = folder / "q_phi_response.csv"
    correlations = {"lensed": [], "unlensed": []}
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("ell_min", "ell_max", "phi_phi_cl", "lensed_q_phi_cl",
                         "unlensed_q_phi_cl", "lensed_q_phi_correlation",
                         "unlensed_q_phi_correlation"))
        for band, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            a, b = results["lensed"][0][band], results["unlensed"][0][band]
            ra = a / np.sqrt(results["lensed"][1][band] * phi_binned[band])
            rb = b / np.sqrt(results["unlensed"][1][band] * phi_binned[band])
            correlations["lensed"].append(ra)
            correlations["unlensed"].append(rb)
            writer.writerow((lo, hi - 1, phi_binned[band], a, b, ra, rb))
            print(f"{lo}-{hi - 1}: r(q_lensed,phi)={ra:+.4f}, r(q_unlensed,phi)={rb:+.4f}")
    centers = (edges[:-1] + edges[1:] - 1) / 2
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for kind, marker, color in (("lensed", "o", "#bb4b2f"),
                                ("unlensed", "s", "#24718c")):
        ax.plot(centers, correlations[kind], marker=marker, color=color,
                lw=1.4, ms=6, label=kind)
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xticks(centers, [f"{lo}–{hi - 1}" for lo, hi in zip(edges[:-1], edges[1:])])
    ax.set_xlabel(r"Output multipole $\ell$ bin")
    ax.set_ylabel(r"Correlation of $q$ with input $\phi$")
    ax.set_title("One ACT-like paired realization", loc="left")
    ax.legend(frameon=False)
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(folder / "q_phi_response.png", dpi=190)
    plt.close(fig)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
