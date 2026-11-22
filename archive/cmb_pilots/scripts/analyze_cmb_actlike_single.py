#!/usr/bin/env python3
"""Plot one ACT-like lensed/unlensed CMB pair in raw bandpower units."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("archive/cmb_pilots/run_cmb_actlike_single_2048"))
    args = parser.parse_args()
    meta = json.loads((args.input_dir / "benchmark.json").read_text())
    if meta["nsims"] != 1:
        parser.error("expected exactly one paired seed")
    rows = list(csv.DictReader((args.input_dir / "bandpowers.csv").open(newline="")))
    if len(rows) != 2 * meta["nbins"]:
        raise ValueError("missing or duplicate bin rows")
    by_kind = {
        kind: sorted((row for row in rows if row["kind"] == kind), key=lambda r: int(r["bin"]))
        for kind in ("lensed", "unlensed")
    }
    if any(len(v) != meta["nbins"] for v in by_kind.values()):
        raise ValueError("each kind must contain every bin")
    bounds = np.array([(int(r["ell_min"]), int(r["ell_max"])) for r in by_kind["lensed"]])
    centers = bounds.mean(axis=1)
    fields = ("marked_minus_phase_controls", "connected_q_minus_phase_controls",
              "connected_q_component_minus_phase_controls")
    residuals = {
        kind: np.array([[float(row[field]) for field in fields] for row in data])
        for kind, data in by_kind.items()
    }
    delta = residuals["lensed"] - residuals["unlensed"]
    scale = 1e7
    fig, axes = plt.subplots(3, 1, figsize=(8.1, 8.8), sharex=True)
    colors = {"lensed": "#bb4b2f", "unlensed": "#24718c"}
    for index, title in enumerate((r"Marked: $\Delta C_\ell^{\delta_m\delta_m}$",
                                   r"Quadratic: $\Delta(C_\ell^{qq}-W_\ell)$")):
        ax = axes[index]
        for kind, marker in (("lensed", "o"), ("unlensed", "s")):
            ax.plot(centers, scale * residuals[kind][:, index], marker=marker,
                    color=colors[kind], lw=1.3, ms=5, label=kind)
        ax.axhline(0, color="0.35", lw=0.8)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_ylabel(r"Residual $C_\ell\times10^7$")
        ax.grid(alpha=0.18)
    axes[0].legend(frameon=False, ncol=2)
    common_lim = 1.12 * max(
        np.max(np.abs(scale * residuals[kind][:, index]))
        for kind in residuals for index in (0, 1)
    )
    for ax in axes[:2]:
        ax.set_ylim(-common_lim, common_lim)
    for label, values, marker, color in (
        ("Marked", delta[:, 0], "o", "#202c40"),
        (r"Normalized $q$ connected part", delta[:, 2], "s", "#9d659f"),
        ("Remainder", delta[:, 0] - delta[:, 2], "^", "#42a07a"),
    ):
        axes[2].plot(centers, scale * values, marker=marker, color=color,
                     lw=1.3, ms=5, label=label)
    axes[2].axhline(0, color="0.35", lw=0.8)
    axes[2].set_title("Lensed minus unlensed, after each phase control", loc="left", fontsize=11)
    axes[2].set_ylabel(r"Paired difference $\times10^7$")
    axes[2].set_xlabel(r"Output multipole $\ell$ (bin centers)")
    axes[2].legend(frameon=False, ncol=3, fontsize=8)
    axes[2].grid(alpha=0.18)
    axes[2].set_xticks(centers, [f"{lo}–{hi}" for lo, hi in bounds])
    fig.suptitle("ACT-like CMB mark: one paired realization, Nside=2048", fontsize=13)
    fig.text(0.5, 0.01, "Temperature modes 600–3000; Gaussian smoothing σ=2′. No ensemble uncertainty.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    output = args.input_dir / "single_pair_results.png"
    fig.savefig(output, dpi=190)
    plt.close(fig)
    table = args.input_dir / "single_pair_contrasts.csv"
    with table.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("ell_min", "ell_max", "marked_pair_difference",
                         "connected_q_pair_difference", "normalized_q_connected_pair_difference",
                         "marked_remainder"))
        for (lo, hi), values in zip(bounds, delta):
            writer.writerow((lo, hi, values[0], values[1], values[2], values[0] - values[2]))
    print(f"Wrote {output} and {table}")


if __name__ == "__main__":
    main()
