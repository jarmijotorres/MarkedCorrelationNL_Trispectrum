#!/usr/bin/env python3
"""Benchmark paired CMB lensing/null spectra for the exact p=2 density mark.

Each seed produces a lensed map, its exact paired unlensed primary, and
phase-randomized controls. Full maps are discarded after binned descriptors
are recorded. The input spectra can be toy spectra or a Planck-like CAMB model.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import healpy as hp
import numpy as np
from lenspyx.lensing import synfast

from cmb.core import bin_spectrum, phase_randomize_alm, toy_spectra


def legendre_matrix(mu: np.ndarray, lmax: int) -> np.ndarray:
    """Return P_ell(mu) with ell on the first axis."""
    matrix = np.empty((lmax + 1, mu.size))
    matrix[0] = 1
    matrix[1] = mu
    for ell in range(2, lmax + 1):
        matrix[ell] = ((2 * ell - 1) * mu * matrix[ell - 1] - (ell - 1) * matrix[ell - 2]) / ell
    return matrix


def fixed_bin_edges(lmax: int, nbins: int) -> np.ndarray:
    """Nonempty, approximately equal-width integer multipole bins."""
    if nbins > lmax - 1:
        raise ValueError("nbins cannot exceed the number of multipoles from 2 through lmax")
    return np.rint(np.linspace(2, lmax + 1, nbins + 1)).astype(int)


def descriptors(
    map_values: np.ndarray, nside: int, lmax: int, window: np.ndarray,
    quadrature: tuple[np.ndarray, np.ndarray],
    normalize: bool = True,
    phi_alm: np.ndarray | None = None,
    raw_spectra: dict[str, np.ndarray] | None = None,
    output_lmax: int | None = None,
    wick_spectrum: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float, np.ndarray | None]:
    """Compute the exact sign-paired p=2 marked spectrum and Q projection.

    With F_+=(1+t)(1+s)^2 and F_-=(1-t)(1-s)^2, write
    F_±=1+y±x, x=t+2s+t*s², y=2ts+s². Their distinct means μ_± must
    each divide the corresponding centered field. The coefficient of C_qq
    in the sign-averaged spectrum is 2(μ_+^-2+μ_-^-2), q=ts-<ts>.
    """
    # Internal t is the toy temperature contrast δ_T, not temperature in K.
    t = np.asarray(map_values, dtype=float)
    t = t - t.mean()
    if normalize:
        t = t / t.std()
    t_alm = hp.map2alm(t, lmax=lmax, iter=3)
    s = hp.alm2map(hp.almxfl(t_alm, window), nside, lmax=lmax)
    ts = t * s
    q = ts - ts.mean()
    s2 = s * s
    x = t + 2 * s + t * s2
    y = 2 * ts + s2
    mu_plus = 1 + x.mean() + y.mean()
    mu_minus = 1 - x.mean() + y.mean()
    if abs(mu_plus) < 1e-12 or abs(mu_minus) < 1e-12:
        raise ValueError("The mean of a marked field is too close to zero")
    x -= x.mean()
    y -= y.mean()
    if output_lmax is None:
        output_lmax = lmax
    q_alm = hp.map2alm(q, lmax=output_lmax, iter=3)
    x_alm = hp.map2alm(x, lmax=output_lmax, iter=3)
    y_alm = hp.map2alm(y, lmax=output_lmax, iter=3)
    cl_tt = hp.alm2cl(t_alm)
    cl_qq = hp.alm2cl(q_alm)
    if raw_spectra is not None:
        raw_spectra.update(tt=cl_tt, qq=cl_qq)
    sum_inverse_means = 1 / mu_plus**2 + 1 / mu_minus**2
    marked = (0.5 * sum_inverse_means * (hp.alm2cl(x_alm) + hp.alm2cl(y_alm))
              + (1 / mu_plus**2 - 1 / mu_minus**2) * hp.alm2cl(x_alm, y_alm))
    ell = np.arange(lmax + 1)
    xi_basis, integration_basis = quadrature
    if wick_spectrum is None:
        xi_tt = cl_tt @ xi_basis
        xi_ss = (cl_tt * window**2) @ xi_basis
        xi_ts = (cl_tt * window) @ xi_basis
        wick = integration_basis @ (xi_tt * xi_ss + xi_ts**2)
    else:
        wick = wick_spectrum
    connected_q = cl_qq[:integration_basis.shape[0]] - wick
    q_coefficient = 2 * sum_inverse_means
    connected_q_component = q_coefficient * connected_q
    if phi_alm is not None:
        ell_out, emm_out = hp.Alm.getlm(output_lmax)
        phi_out = phi_alm[hp.Alm.getidx(lmax, ell_out, emm_out)]
        q_phi = hp.alm2cl(q_alm, phi_out)
    else:
        q_phi = None
    return (marked, connected_q, connected_q_component, cl_tt, t_alm,
            float(mu_plus), float(mu_minus), float(q_coefficient), q_phi)


def main() -> None:
    process_start = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside", type=int, default=128)
    parser.add_argument("--lmax", type=int, default=128)
    parser.add_argument("--nsims", type=int, default=10)
    parser.add_argument("--surrogates", type=int, default=2)
    parser.add_argument("--nbins", type=int, default=12)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--smoothing-degrees", type=float, default=2 / 60,
        help="Gaussian sigma of the smoothed temperature contrast in degrees (default: 2 arcmin)",
    )
    parser.add_argument("--temperature-lmin", type=int, default=2)
    parser.add_argument("--spectra", choices=("toy", "camb"), default="toy")
    parser.add_argument("--bin-edges", type=int, nargs="+")
    parser.add_argument("--save-unbinned", action="store_true")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output-dir", type=Path, default=Path("run_cmb_benchmark_128_p2"))
    args = parser.parse_args()
    if args.nside < 1 or args.lmax < 2 or args.lmax > 3 * args.nside - 1:
        parser.error("require nside >= 1 and 2 <= lmax <= 3*nside-1")
    if args.nsims < 1 or args.surrogates < 1 or args.threads < 1:
        parser.error("nsims, surrogates and threads must be positive")
    if not (2 <= args.temperature_lmin <= args.lmax):
        parser.error("temperature-lmin must lie between 2 and lmax")
    if args.bin_edges is None:
        edges = fixed_bin_edges(args.lmax, args.nbins)
    else:
        edges = np.asarray(args.bin_edges)
        if (len(edges) < 2 or edges[0] < 2 or edges[-1] > args.lmax + 1
                or np.any(np.diff(edges) <= 0)):
            parser.error("bin-edges must increase from at least 2 to at most lmax+1")
        args.nbins = len(edges) - 1
    if args.spectra == "toy":
        ell, tt, pp = toy_spectra(args.lmax)
    else:
        import camb

        parameters = camb.CAMBparams()
        parameters.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.12, mnu=0.06, tau=0.054)
        parameters.InitPower.set_params(As=2.1e-9, ns=0.965)
        parameters.set_for_lmax(args.lmax, lens_potential_accuracy=1)
        results = camb.get_results(parameters)
        tt = results.get_cmb_power_spectra(
            parameters, CMB_unit="muK", raw_cl=True,
        )["unlensed_scalar"][:args.lmax + 1, 0]
        pp = results.get_lens_potential_cls(lmax=args.lmax, raw_cl=True)[:, 0]
        ell = np.arange(args.lmax + 1)
    window = np.exp(-0.5 * ell * (ell + 1) * np.deg2rad(args.smoothing_degrees) ** 2)
    output_lmax = int(edges[-1]) - 1
    from scipy.special import roots_legendre

    mu, weights = roots_legendre(args.lmax + (output_lmax + 1) // 2 + 1)
    polynomials = legendre_matrix(mu, args.lmax)
    xi_basis = (2 * ell[:, None] + 1) * polynomials / (4 * np.pi)
    integration_basis = 2 * np.pi * polynomials[:output_lmax + 1] * weights[None, :]
    quadrature = (xi_basis, integration_basis)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "input_spectra.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("ell", "unlensed_tt_cl_uK2", "phi_phi_cl"))
        writer.writerows(zip(ell, tt, pp))
    rng = np.random.default_rng(args.seed)
    rows = []
    all_times = []
    max_phase_power_mismatch = 0.0
    unbinned = {}
    start = time.perf_counter()
    setup_seconds = start - process_start
    stage_times = {"generation": 0.0, "filtering": 0.0,
                   "descriptors": 0.0, "controls": 0.0}
    for index in range(args.nsims):
        seed = int(rng.integers(0, 2**31 - 1))
        tick = time.perf_counter()
        maps, (alms, labels) = synfast(
            {"tt": tt, "pp": pp}, lmax=args.lmax,
            geometry=("healpix", {"nside": args.nside}),
            seed=seed, alm=True, nthreads=args.threads, verbose=False,
        )
        phi_power = hp.alm2cl(alms[labels.index("p")])
        if args.save_unbinned:
            unbinned[f"pair_{index}_phi_phi"] = phi_power
        unlensed = hp.alm2map(alms[labels.index("t")], args.nside, lmax=args.lmax)
        # Common random phases make the lensed/unlensed controls a paired
        # comparison, reducing avoidable surrogate Monte Carlo variance.
        phase_seeds = rng.integers(0, 2**31 - 1, size=args.surrogates)
        generated = time.perf_counter()
        stage_times["generation"] += generated - tick
        for kind, values in (("lensed", maps["T"]), ("unlensed", unlensed)):
            stage_tick = time.perf_counter()
            if args.temperature_lmin > 2:
                filtered_alm = hp.map2alm(values, lmax=args.lmax, iter=1)
                filtered_alm = hp.almxfl(
                    filtered_alm, np.where(ell >= args.temperature_lmin, 1.0, 0.0),
                )
                values = hp.alm2map(filtered_alm, args.nside, lmax=args.lmax)
                del filtered_alm
            stage_times["filtering"] += time.perf_counter() - stage_tick
            raw_spectra = {} if args.save_unbinned else None
            stage_tick = time.perf_counter()
            marked, trispectrum, q_component, power, field_alm, mu_plus, mu_minus, q_coefficient, q_phi = descriptors(
                values, args.nside, args.lmax, window, quadrature,
                phi_alm=alms[labels.index("p")],
                raw_spectra=raw_spectra,
                output_lmax=output_lmax,
            )
            stage_times["descriptors"] += time.perf_counter() - stage_tick
            prefix = f"pair_{index}_{kind}"
            if raw_spectra is not None:
                for name, spectrum in raw_spectra.items():
                    unbinned[f"{prefix}_{name}"] = spectrum
                unbinned[f"{prefix}_q_coefficient"] = q_coefficient
            surrogate_marked = []
            surrogate_trispectrum = []
            surrogate_q_component = []
            # Phase randomization retains each alm amplitude, so its Wick
            # spectrum is the same as the original map up to map/alm roundoff.
            original_wick = (raw_spectra["qq"][:output_lmax + 1] - trispectrum
                             if raw_spectra is not None else None)
            stage_tick = time.perf_counter()
            for control_index, phase_seed in enumerate(phase_seeds):
                randomized_alm = phase_randomize_alm(field_alm, args.lmax, phase_seed)
                randomized = hp.alm2map(randomized_alm, args.nside, lmax=args.lmax)
                control_spectra = {} if args.save_unbinned else None
                surrogate = descriptors(
                    randomized, args.nside, args.lmax, window, quadrature,
                    normalize=False, raw_spectra=control_spectra,
                    output_lmax=output_lmax, wick_spectrum=original_wick,
                )
                if control_spectra is not None:
                    control_prefix = f"{prefix}_control_{control_index}"
                    for name, spectrum in control_spectra.items():
                        unbinned[f"{control_prefix}_{name}"] = spectrum
                    unbinned[f"{control_prefix}_q_coefficient"] = surrogate[7]
                    unbinned[f"{control_prefix}_marked"] = surrogate[0]
                max_phase_power_mismatch = max(
                    max_phase_power_mismatch,
                    float(np.max(
                        np.abs(power[args.temperature_lmin:] - surrogate[3][args.temperature_lmin:])
                        / np.maximum(power[args.temperature_lmin:],
                                     1e-12 * np.max(power[args.temperature_lmin:])),
                    )),
                )
                surrogate_marked.append(surrogate[0])
                surrogate_trispectrum.append(surrogate[1])
                surrogate_q_component.append(surrogate[2])
            stage_times["controls"] += time.perf_counter() - stage_tick
            marked_residual = marked - np.mean(surrogate_marked, axis=0)
            trispectrum_residual = trispectrum - np.mean(surrogate_trispectrum, axis=0)
            q_component_residual = q_component - np.mean(surrogate_q_component, axis=0)
            if args.save_unbinned:
                unbinned[f"pair_{index}_{kind}_marked"] = marked
                unbinned[f"pair_{index}_{kind}_connected_q"] = trispectrum
                unbinned[f"pair_{index}_{kind}_q_component"] = q_component
                unbinned[f"pair_{index}_{kind}_marked_residual"] = marked_residual
                unbinned[f"pair_{index}_{kind}_connected_q_residual"] = trispectrum_residual
                unbinned[f"pair_{index}_{kind}_q_component_residual"] = q_component_residual
                unbinned[f"pair_{index}_{kind}_q_phi"] = q_phi
            bandpowers = [
                bin_spectrum(spectrum, edges)
                for spectrum in (marked, trispectrum, q_component, power,
                                 marked_residual, trispectrum_residual, q_component_residual,
                                 q_phi, phi_power)
            ]
            for band, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                binned = [values[band] for values in bandpowers]
                rows.append((index, seed, kind, band, lo, hi - 1, *binned,
                             mu_plus, mu_minus, q_coefficient))
        all_times.append(time.perf_counter() - tick)
        print(f"pair {index + 1}/{args.nsims}: generation {generated - tick:.2f}s, total {all_times[-1]:.2f}s", flush=True)
    headings = (
        "pair", "seed", "kind", "bin", "ell_min", "ell_max",
        "marked_cl", "connected_q_cl", "connected_q_component_cl", "temperature_cl",
        "marked_minus_phase_controls", "connected_q_minus_phase_controls",
        "connected_q_component_minus_phase_controls",
        "q_phi_cl", "phi_phi_cl",
        "mu_plus", "mu_minus", "q_connected_coefficient",
    )
    with (args.output_dir / "bandpowers.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headings)
        writer.writerows(rows)
    if args.save_unbinned:
        np.savez_compressed(args.output_dir / "unbinned_spectra.npz", **unbinned)
    metadata = vars(args).copy()
    metadata["mark_exponent"] = 2
    metadata["smoothing_sigma_arcmin"] = 60 * args.smoothing_degrees
    metadata["smoothing_fwhm_arcmin"] = (2 * np.sqrt(2 * np.log(2))
                                            * 60 * args.smoothing_degrees)
    metadata["smoothing_field"] = "filtered, unit-variance temperature contrast; not phi"
    metadata["smoothing_window"] = "exp[-ell(ell+1)*sigma_radians^2/2]"
    metadata["temperature_contrast"] = "delta_T = (T-mean(T))/std(T) for each map"
    metadata["marked_field"] = "(1+delta_T)(1+s)^2/<(1+delta_T)(1+s)^2>-1"
    metadata["paired_phase_control_seeds"] = True
    if args.spectra == "camb":
        metadata["camb_parameters"] = {
            "H0": 67.4, "ombh2": 0.0224, "omch2": 0.12, "mnu": 0.06,
            "tau": 0.054, "As": 2.1e-9, "ns": 0.965,
            "lens_potential_accuracy": 1,
        }
    metadata["output_dir"] = str(args.output_dir)
    metadata["bin_edges"] = edges.tolist()
    metadata["mean_seconds_per_pair"] = float(np.mean(all_times))
    metadata["total_seconds"] = time.perf_counter() - start
    metadata["setup_seconds"] = setup_seconds
    metadata["stage_seconds"] = stage_times
    metadata["max_phase_power_mismatch"] = max_phase_power_mismatch
    (args.output_dir / "benchmark.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Mean {np.mean(all_times):.2f}s per paired seed; {metadata['total_seconds']:.1f}s total")
    print(f"Maximum realized phase-control C_ell mismatch: {max_phase_power_mismatch:.2e}")
    print(f"Projected serial time: 100 pairs {np.mean(all_times) * 100 / 60:.1f} min; "
          f"1000 pairs {np.mean(all_times) * 1000 / 3600:.1f} h")


if __name__ == "__main__":
    main()
