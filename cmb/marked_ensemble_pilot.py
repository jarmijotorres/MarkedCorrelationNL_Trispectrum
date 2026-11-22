#!/usr/bin/env python3
"""Run and monitor a resumable paired CMB marked-spectrum ensemble."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from cmb.core import bin_spectrum


def load_config(path: Path, count: int | None = None,
                output_dir: Path | None = None) -> dict:
    config = json.loads(path.read_text())
    if "seed_base" in config:
        requested = count if count is not None else config.get("count")
        if requested is None or requested < 1:
            raise ValueError("provide a positive --count or config count")
        config["seeds"] = np.random.default_rng(config["seed_base"]).integers(
            0, 2**31 - 1, size=requested).tolist()
        config["count"] = requested
    elif count is not None:
        if count < 1 or count > len(config["seeds"]):
            raise ValueError("--count must lie within the explicit seed list")
        config["seeds"] = config["seeds"][:count]
    if len(config["seeds"]) != len(set(config["seeds"])) or not config["seeds"]:
        raise ValueError("seeds must be a nonempty list of distinct integers")
    edges = config["bin_edges"]
    if len(edges) < 2 or edges[0] < 2 or edges[-1] > config["lmax"] + 1 or any(
        hi <= lo for lo, hi in zip(edges[:-1], edges[1:])
    ):
        raise ValueError("invalid bin_edges")
    for low, high in config["summary_bands"]:
        if not (edges[0] <= low < high <= edges[-1]):
            raise ValueError("summary band outside output multipoles")
    config["output_dir"] = str((output_dir if output_dir is not None
                                else path.parent.parent / config["output_dir"]).resolve())
    return config


def pair_dir(config: dict, index: int) -> Path:
    return Path(config["output_dir"]) / f"pair_{index:02d}"


def complete(config: dict, index: int) -> bool:
    folder = pair_dir(config, index)
    required = ("benchmark.json", "bandpowers.csv", "unbinned_spectra.npz")
    if not all((folder / name).exists() for name in required):
        return False
    meta = json.loads((folder / "benchmark.json").read_text())
    for key in ("nside", "lmax", "temperature_lmin", "spectra", "surrogates", "threads", "bin_edges"):
        if meta[key] != config[key]:
            raise ValueError(f"{folder}: {key} differs from config")
    if (meta["seed"] != config["seeds"][index] or meta["nsims"] != 1
            or not meta["save_unbinned"]
            or not np.isclose(meta["smoothing_degrees"], config["smoothing_degrees"])):
        raise ValueError(f"{folder}: run settings differ from config")
    return True


def band(cl: np.ndarray, low: int, high: int) -> float:
    return float(bin_spectrum(cl, np.array([low, high]))[0])


def write_progress(config: dict, records: list[dict]) -> None:
    root = Path(config["output_dir"])
    if records:
        temporary = root / "summary.csv.tmp"
        with temporary.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        temporary.replace(root / "summary.csv")
    progress = {"completed": len(records), "requested": len(config["seeds"]),
                "summary_bands": {}}
    for low, high in config["summary_bands"]:
        label = f"L{low}_{high}"
        values = np.array([row[f"marked_{label}"] for row in records])
        progress["summary_bands"][label] = {
            "mean": float(values.mean()) if len(values) else None,
            "sample_sd": float(values.std(ddof=1)) if len(values) > 1 else None,
            "positive": int(np.count_nonzero(values > 0)),
        }
    temporary = root / "progress.json.tmp"
    temporary.write_text(json.dumps(progress, indent=2) + "\n")
    temporary.replace(root / "progress.json")


def pair_record(config: dict, index: int) -> dict:
    with np.load(pair_dir(config, index) / "unbinned_spectra.npz") as saved:
        paired = (saved["pair_0_lensed_marked_residual"]
                  - saved["pair_0_unlensed_marked_residual"])
        controls = np.array([
            saved[f"pair_0_lensed_control_{j}_marked"]
            - saved[f"pair_0_unlensed_control_{j}_marked"]
            for j in range(config["surrogates"])
        ])
        raw = saved["pair_0_lensed_marked"] - saved["pair_0_unlensed_marked"]
        if not np.allclose(paired, raw - controls.mean(axis=0),
                           rtol=1e-9, atol=1e-13):
            raise ValueError(f"pair {index}: phase-control reconstruction failed")
        q_paired = (saved["pair_0_lensed_q_component_residual"]
                    - saved["pair_0_unlensed_q_component_residual"])
        record = {"pair": index, "master_seed": config["seeds"][index]}
        for low, high in config["summary_bands"]:
            label = f"L{low}_{high}"
            record[f"marked_{label}"] = band(paired, low, high)
            record[f"q_component_{label}"] = band(q_paired, low, high)
            control_bands = [band(c, low, high) for c in controls]
            record[f"control_contrast_sd_{label}"] = float(
                np.std(control_bands, ddof=1)
            ) if len(control_bands) > 1 else float("nan")
        return record


def diagnose(config: dict, make_plot: bool = True) -> list[dict]:
    root = Path(config["output_dir"])
    root.mkdir(parents=True, exist_ok=True)
    records = []
    spectra = []
    lensed_spectra = []
    unlensed_spectra = []
    for index, master_seed in enumerate(config["seeds"]):
        if not complete(config, index):
            continue
        with np.load(pair_dir(config, index) / "unbinned_spectra.npz") as saved:
            lensed = saved["pair_0_lensed_marked_residual"]
            unlensed = saved["pair_0_unlensed_marked_residual"]
            paired = lensed - unlensed
            controls = np.array([
                saved[f"pair_0_lensed_control_{j}_marked"]
                - saved[f"pair_0_unlensed_control_{j}_marked"]
                for j in range(config["surrogates"])
            ])
            # The paired residual must equal paired raw marked minus the
            # mean of the matched phase-control contrasts.
            raw = saved["pair_0_lensed_marked"] - saved["pair_0_unlensed_marked"]
            if not np.allclose(paired, raw - controls.mean(axis=0), rtol=1e-9, atol=1e-13):
                raise ValueError(f"pair {index}: phase-control reconstruction failed")
            q_paired = (saved["pair_0_lensed_q_component_residual"]
                        - saved["pair_0_unlensed_q_component_residual"])
            record = {"pair": index, "master_seed": master_seed}
            for low, high in config["summary_bands"]:
                label = f"L{low}_{high}"
                record[f"marked_{label}"] = band(paired, low, high)
                record[f"q_component_{label}"] = band(q_paired, low, high)
                control_bands = [band(c, low, high) for c in controls]
                record[f"control_contrast_sd_{label}"] = float(
                    np.std(control_bands, ddof=1)
                ) if len(control_bands) > 1 else float("nan")
            records.append(record)
            spectra.append(paired.copy())
            lensed_spectra.append(lensed.copy())
            unlensed_spectra.append(unlensed.copy())

    write_progress(config, records)
    if records and make_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        edges = np.asarray(config["bin_edges"])
        centers = bin_spectrum(np.arange(config["lmax"] + 1), edges)
        binned = np.array([bin_spectrum(cl, edges) for cl in spectra])
        binned_lensed = np.array([bin_spectrum(cl, edges) for cl in lensed_spectra])
        binned_unlensed = np.array([bin_spectrum(cl, edges) for cl in unlensed_spectra])
        binned_q = []
        for record in records:
            with np.load(pair_dir(config, record["pair"]) / "unbinned_spectra.npz") as saved:
                q = (saved["pair_0_lensed_q_component_residual"]
                     - saved["pair_0_unlensed_q_component_residual"])
                binned_q.append(bin_spectrum(q, edges))
        binned_q = np.asarray(binned_q)
        np.savez_compressed(
            root / "ensemble_vectors.npz", bin_edges=edges, centers=centers,
            paired_marked=binned, paired_q_component=binned_q,
            paired_remainder=binned - binned_q,
            lensed=binned_lensed, unlensed=binned_unlensed,
            sample_covariance=np.cov(binned, rowvar=False) if len(records) > 1
            else np.full((len(centers), len(centers)), np.nan),
            master_seeds=np.asarray([r["master_seed"] for r in records]),
        )
        fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                 constrained_layout=True)
        if len(records) <= 30:
            for values in binned_lensed:
                axes[0].plot(centers, values * 1e7, lw=0.9, alpha=0.32,
                             color="#bb4b2f")
            for values in binned_unlensed:
                axes[0].plot(centers, values * 1e7, lw=0.9, alpha=0.32,
                             color="#24718c")
        else:
            for values, color in ((binned_lensed, "#bb4b2f"),
                                  (binned_unlensed, "#24718c")):
                lower, upper = np.quantile(values, [0.16, 0.84], axis=0)
                axes[0].fill_between(centers, lower * 1e7, upper * 1e7,
                                     alpha=0.15, color=color)
        axes[0].plot(centers, binned_lensed.mean(axis=0) * 1e7, lw=2.4,
                     color="#bb4b2f", label="lensed mean")
        axes[0].plot(centers, binned_unlensed.mean(axis=0) * 1e7, lw=2.4,
                     color="#24718c", label="unlensed mean")
        axes[0].set(ylabel=r"$10^7 C_L^{\rm mark,res}$",
                    title=("Phase-control residuals: individual skies"
                           if len(records) <= 30 else
                           "Phase-control residuals: 16–84% sky range"))
        unlensed_spread = binned_unlensed.std(axis=0, ddof=1) if len(records) > 1 else np.zeros_like(centers)
        axes[1].fill_between(centers, -unlensed_spread * 1e7,
                             unlensed_spread * 1e7, color="#24718c",
                             alpha=0.2, label=r"unlensed spread ($\pm 1$ sample SD)")
        axes[1].plot(centers, binned.mean(axis=0) * 1e7, lw=2.5,
                     color="#202c40", label="paired mean")
        axes[1].set(xlabel="Output multipole L",
                    ylabel=r"$10^7\Delta C_L^{\rm mark,res}$",
                    title="Paired residual: lensed minus unlensed")
        for ax in axes:
            ax.axhline(0, color="0.4", lw=0.8)
            ax.grid(alpha=0.2)
            ax.legend(frameon=False, fontsize=9)
        fig.suptitle(f"{len(records)}/{len(config['seeds'])} paired skies; "
                     f"{config['surrogates']} phase controls per map")
        fig.savefig(root / "marked_diagnostic.png", dpi=180)
        plt.close(fig)
    return records


def run(config: dict) -> None:
    root = Path(config["output_dir"])
    root.mkdir(parents=True, exist_ok=True)
    (root / "config_snapshot.json").write_text(json.dumps(config, indent=2) + "\n")
    records = diagnose(config, make_plot=False)
    for index, seed in enumerate(config["seeds"]):
        if complete(config, index):
            print(f"pair {index + 1}/{len(config['seeds'])}: already complete", flush=True)
            continue
        folder = pair_dir(config, index)
        folder.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable, "-m", "cmb.benchmark_cmb_ensemble",
            "--nside", str(config["nside"]), "--lmax", str(config["lmax"]),
            "--nsims", "1", "--surrogates", str(config["surrogates"]),
            "--seed", str(seed), "--smoothing-degrees", str(config["smoothing_degrees"]),
            "--temperature-lmin", str(config["temperature_lmin"]),
            "--spectra", config["spectra"], "--threads", str(config["threads"]),
            "--bin-edges", *map(str, config["bin_edges"]), "--save-unbinned",
            "--output-dir", str(folder),
        ]
        print(f"pair {index + 1}/{len(config['seeds'])}: starting master seed {seed}", flush=True)
        env = os.environ.copy()
        env.setdefault("OPENBLAS_NUM_THREADS", "1")
        env.setdefault("OMP_NUM_THREADS", str(config["threads"]))
        with (folder / "run.log").open("w") as log:
            process = subprocess.Popen(command, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True, env=env)
            assert process.stdout is not None
            for line in process.stdout:
                print(line, end="", flush=True)
                log.write(line)
                log.flush()
            code = process.wait()
        if code:
            raise SystemExit(f"pair {index + 1} failed with exit code {code}; see {folder / 'run.log'}")
        if not complete(config, index):
            raise RuntimeError(f"pair {index + 1} exited without complete outputs")
        records.append(pair_record(config, index))
        records.sort(key=lambda item: item["pair"])
        write_progress(config, records)
        low, high = config["summary_bands"][0]
        print(f"completed {len(records)}/{len(config['seeds'])}; "
              f"paired marked C_L[{low},{high})={records[-1][f'marked_L{low}_{high}']:.3e}",
              flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "status", "diagnose"))
    parser.add_argument("--config", type=Path, default=Path("configs/cmb_marked_10.json"))
    parser.add_argument("--count", type=int,
                        help="number of pairs; generated seeds keep the same prefix when extended")
    parser.add_argument("--output-dir", type=Path,
                        help="override output directory, e.g. a cluster scratch path")
    args = parser.parse_args()
    config = load_config(args.config.resolve(), args.count, args.output_dir)
    if args.action == "run":
        run(config)
    elif args.action == "status":
        root = Path(config["output_dir"])
        done = sum(complete(config, i) for i in range(len(config["seeds"])))
        print(f"{done}/{len(config['seeds'])} paired skies complete")
        progress_file = root / "progress.json"
        if progress_file.exists():
            progress = json.loads(progress_file.read_text())
            if progress["completed"] == done:
                for label, result in progress["summary_bands"].items():
                    if result["mean"] is not None:
                        print(f"{label}: mean {result['mean']:.3e}, "
                              f"positive {result['positive']}/{done}")
    else:
        records = diagnose(config)
        print(f"{len(records)}/{len(config['seeds'])} paired skies complete")
        for low, high in config["summary_bands"]:
            values = [row[f"marked_L{low}_{high}"] for row in records]
            if values:
                print(f"{low} ≤ L < {high}: mean {np.mean(values):.3e}, "
                      f"positive {sum(v > 0 for v in values)}/{len(values)}")


if __name__ == "__main__":
    main()
