"""Replay saved CMB pair seeds to measure linear-mark auto and Ebina cross terms.

The main run discarded maps. This follow-up deterministically recreates each
map and its ten phase controls from the saved configuration, then writes only
the new binned descriptors. It never changes the original pair products.

python -m cmb.linear_mark_followup run --config configs/cmb_marked_log8_1000.json --workers 64
python -m cmb.linear_mark_followup diagnose --config configs/cmb_marked_log8_1000.json
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from cmb.core import bin_spectrum, phase_randomize_alm
from cmb.marked_ensemble_pilot import complete, load_config as load_ensemble_config, pair_dir


def load_config(path):
    """Load the ensemble configuration and retain its path for child workers."""
    path = Path(path).resolve()
    config = load_ensemble_config(path)
    config["config_path"] = str(path)
    return config


def spectra_for_map(values, nside, lmax, output_lmax, window, normalize):
    """Return sign-paired linear auto, unpaired Ebina cubic, and temperature alm.

    Ebina is the +mark unmarked--marked cross after subtracting the exact
    two-point contribution: C_L[t, ts-<ts>]/< (1+t)(1+s) >.
    The sign-paired auto averages the + and - linear marks, each normalized
    by its own realized mean.
    """
    import healpy as hp

    t = np.asarray(values, dtype=float)
    t = t - t.mean()
    if normalize:
        t = t / t.std()
    t_alm = hp.map2alm(t, lmax=lmax, iter=3)
    s = hp.alm2map(hp.almxfl(t_alm, window), nside, lmax=lmax)
    ts = t * s
    u = ts.mean()
    x = t + s
    x_mean = x.mean()
    x = x - x_mean
    q = ts - u
    mu_plus = 1 + x_mean + u
    mu_minus = 1 - x_mean + u
    if min(abs(mu_plus), abs(mu_minus)) < 1e-12:
        raise ValueError("Linear-mark normalization is too close to zero")
    x_alm = hp.map2alm(x, lmax=output_lmax, iter=3)
    q_alm = hp.map2alm(q, lmax=output_lmax, iter=3)
    c_xx = hp.alm2cl(x_alm)
    c_qq = hp.alm2cl(q_alm)
    c_xq = hp.alm2cl(x_alm, q_alm)
    auto = 0.5 * ((c_xx + c_qq + 2 * c_xq) / mu_plus**2
                  + (c_xx + c_qq - 2 * c_xq) / mu_minus**2)
    ell_out, emm_out = hp.Alm.getlm(output_lmax)
    t_out = t_alm[hp.Alm.getidx(lmax, ell_out, emm_out)]
    ebina = hp.alm2cl(t_out, q_alm) / mu_plus
    return auto, ebina, t_alm


def replay_pair(config, index):
    import healpy as hp
    from lenspyx.lensing import synfast

    folder = pair_dir(config, index)
    if not complete(config, index):
        raise ValueError(f"Original pair {index} is incomplete")
    # This CSV is the exact input TT/phi spectrum used by the original run.
    spectra = np.loadtxt(folder / "input_spectra.csv", delimiter=",", skiprows=1)
    lmax = config["lmax"]
    tt, pp = spectra[:lmax + 1, 1], spectra[:lmax + 1, 2]
    rng = np.random.default_rng(config["seeds"][index])
    synfast_seed = int(rng.integers(0, 2**31 - 1))
    maps, (alms, labels) = synfast(
        {"tt": tt, "pp": pp}, lmax=lmax,
        geometry=("healpix", {"nside": config["nside"]}),
        seed=synfast_seed, alm=True, nthreads=config["threads"], verbose=False,
    )
    unlensed = hp.alm2map(alms[labels.index("t")], config["nside"], lmax=lmax)
    phase_seeds = rng.integers(0, 2**31 - 1, size=config["surrogates"])
    ell = np.arange(lmax + 1)
    window = np.exp(-0.5 * ell * (ell + 1)
                    * np.deg2rad(config["smoothing_degrees"]) ** 2)
    passband = np.where(ell >= config["temperature_lmin"], 1.0, 0.0)
    edges = np.asarray(config["bin_edges"])
    output_lmax = int(edges[-1]) - 1
    output = {}
    with np.load(folder / "unbinned_spectra.npz") as saved:
        for kind, values in (("lensed", maps["T"]), ("unlensed", unlensed)):
            if config["temperature_lmin"] > 2:
                filtered_alm = hp.map2alm(values, lmax=lmax, iter=1)
                values = hp.alm2map(hp.almxfl(filtered_alm, passband),
                                    config["nside"], lmax=lmax)
            auto, ebina, t_alm = spectra_for_map(
                values, config["nside"], lmax, output_lmax, window, normalize=True
            )
            # Detect a replay mismatch before trusting the new descriptors.
            reference_tt = saved[f"pair_0_{kind}_tt"]
            if not np.allclose(hp.alm2cl(t_alm), reference_tt,
                               rtol=2e-5, atol=1e-12 * np.max(reference_tt)):
                raise ValueError(f"Pair {index} {kind}: temperature replay differs from saved spectrum")
            controls_auto, controls_ebina = [], []
            for phase_seed in phase_seeds:
                randomized_alm = phase_randomize_alm(t_alm, lmax, int(phase_seed))
                randomized = hp.alm2map(randomized_alm, config["nside"], lmax=lmax)
                control_auto, control_ebina, _ = spectra_for_map(
                    randomized, config["nside"], lmax, output_lmax, window,
                    normalize=False,
                )
                controls_auto.append(bin_spectrum(control_auto, edges))
                controls_ebina.append(bin_spectrum(control_ebina, edges))
            auto = bin_spectrum(auto, edges)
            ebina = bin_spectrum(ebina, edges)
            controls_auto = np.asarray(controls_auto)
            controls_ebina = np.asarray(controls_ebina)
            output[f"{kind}_linear_auto"] = auto
            output[f"{kind}_linear_auto_controls"] = controls_auto
            output[f"{kind}_linear_auto_residual"] = auto - controls_auto.mean(axis=0)
            output[f"{kind}_ebina_bispectrum"] = ebina
            output[f"{kind}_ebina_controls"] = controls_ebina
            output[f"{kind}_ebina_residual"] = ebina - controls_ebina.mean(axis=0)
    output["paired_linear_auto_residual"] = (
        output["lensed_linear_auto_residual"] - output["unlensed_linear_auto_residual"]
    )
    output["paired_ebina_residual"] = (
        output["lensed_ebina_residual"] - output["unlensed_ebina_residual"]
    )
    output["master_seed"] = config["seeds"][index]
    output["bin_edges"] = edges
    destination = folder / "linear_mark_followup.npz"
    temporary = folder / "linear_mark_followup.npz.tmp"
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **output)
    temporary.replace(destination)


def validated_followup(config, index):
    path = pair_dir(config, index) / "linear_mark_followup.npz"
    if not path.exists():
        return False
    with np.load(path) as data:
        return (int(data["master_seed"]) == config["seeds"][index]
                and np.array_equal(data["bin_edges"], config["bin_edges"]))


def run(config, workers):
    pending = [i for i in range(config["count"]) if not validated_followup(config, i)]
    print(f"{config['count'] - len(pending)} complete; {len(pending)} remaining", flush=True)
    def launch(index):
        env = os.environ.copy()
        env.update(OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   NUMEXPR_NUM_THREADS="1", OMP_NUM_THREADS=str(config["threads"]))
        command = [sys.executable, "-m", "cmb.linear_mark_followup", "pair",
                   "--config", config["config_path"], "--index", str(index)]
        folder = pair_dir(config, index)
        with (folder / "linear_mark_followup.log").open("w") as log:
            result = subprocess.run(command, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, check=False)
        if result.returncode:
            raise RuntimeError(f"Pair {index} failed; see {folder / 'linear_mark_followup.log'}")
        return index
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(launch, i) for i in pending]
        for future in as_completed(futures):
            index = future.result()
            print(f"linear follow-up {index + 1}/{config['count']} complete", flush=True)


def diagnose(config):
    fields = ("lensed_linear_auto_residual", "unlensed_linear_auto_residual",
              "paired_linear_auto_residual", "lensed_ebina_bispectrum",
              "unlensed_ebina_bispectrum", "lensed_ebina_residual",
              "unlensed_ebina_residual", "paired_ebina_residual")
    samples = {key: [] for key in fields}
    seeds = []
    for index in range(config["count"]):
        if not validated_followup(config, index):
            continue
        with np.load(pair_dir(config, index) / "linear_mark_followup.npz") as data:
            for key in fields:
                samples[key].append(data[key])
            seeds.append(int(data["master_seed"]))
    if not seeds:
        raise ValueError("No completed linear-mark follow-up pairs")
    root = Path(config["output_dir"])
    output = root / "linear_mark_followup_vectors.npz"
    np.savez_compressed(output, bin_edges=config["bin_edges"], master_seeds=seeds,
                        **{key: np.asarray(value) for key, value in samples.items()})
    print(f"Saved {len(seeds)}/{config['count']} pairs to {output}")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "pair", "diagnose"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--index", type=int)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    config = load_config(args.config.resolve())
    if args.action == "pair":
        if args.index is None or not 0 <= args.index < config["count"]:
            parser.error("pair requires a valid --index")
        replay_pair(config, args.index)
    elif args.action == "run":
        if args.workers < 1:
            parser.error("--workers must be positive")
        run(config, args.workers)
    else:
        diagnose(config)


if __name__ == "__main__":
    main()
