#!/usr/bin/env python3
"""Legacy polynomial-mark comparison on a full-sky CMB temperature map.

This full-sky implementation needs no NaMaster coupling correction.  It
compares a nonlinear marked auto-correlation with a linear-mark (Ebina-style)
cross-correlation, using a harmonic-phase-randomized map with identical C_ell
as the control.  The default input is the lensing-trispectrum toy map created
by ``cmb_lensing_trispectrum_demo.py``. Its polynomial mark differs from the
original p=2 mean-normalized density-style mark; use
``benchmark_cmb_ensemble.py`` for that exact CMB analogue.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def phase_randomize_alm(alm: np.ndarray, lmax: int, seed: int) -> np.ndarray:
    """Randomize HEALPix alm phases while preserving every |a_lm| exactly."""
    import healpy as hp

    rng = np.random.default_rng(seed)
    ell, emm = hp.Alm.getlm(lmax)
    phases = np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, size=alm.size))
    # m=0 coefficients encode real modes, so they may only change sign.
    phases[emm == 0] = rng.choice(np.array([-1.0, 1.0]), size=np.count_nonzero(emm == 0))
    phases[ell == 0] = 1.0
    return alm * phases


def smooth_map(field: np.ndarray, lmax: int, sigma_radians: float) -> np.ndarray:
    """Gaussian-harmonic smoothing of a HEALPix scalar map."""
    import healpy as hp

    alm = hp.map2alm(field, lmax=lmax, iter=3)
    ell = np.arange(lmax + 1)
    window = np.exp(-0.5 * ell * (ell + 1.0) * sigma_radians**2)
    return hp.alm2map(hp.almxfl(alm, window), hp.get_nside(field), lmax=lmax)


def correlation_from_maps(left: np.ndarray, right: np.ndarray, lmax: int, theta: np.ndarray) -> np.ndarray:
    """Return the isotropic angular correlation from the full-sky cross-spectrum."""
    import healpy as hp

    cl = hp.anafast(left, right, lmax=lmax, iter=3)
    legendre_coefficients = (2.0 * np.arange(lmax + 1) + 1.0) * cl / (4.0 * np.pi)
    return np.polynomial.legendre.legval(np.cos(theta), legendre_coefficients)


def descriptors(
    temperature: np.ndarray, lmax: int, smoothing_radians: float, linear_mark: float, quadratic_mark: float,
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return marked, linear-mark, quadratic, and connected-four-point correlations.

    The CMB has no density-pair factor.  We therefore use the generic mark
    m(s)=1+a*s+b*(s^2-<s^2>) with F=t*m(s), and define
    X_E=<t_1 [t*m_linear(s)]_2> for the linear-mark comparison.
    """
    temperature = temperature - temperature.mean()
    smoothed = smooth_map(temperature, lmax, smoothing_radians)
    nonlinear_mark = 1.0 + linear_mark * smoothed + quadratic_mark * (smoothed**2 - np.mean(smoothed**2))
    nonlinear_field = temperature * nonlinear_mark
    nonlinear_field -= nonlinear_field.mean()
    linear_field = temperature * (1.0 + linear_mark * smoothed)
    linear_field -= linear_field.mean()
    quadratic_field = temperature * smoothed
    quadratic_field -= quadratic_field.mean()
    c_quadratic = correlation_from_maps(quadratic_field, quadratic_field, lmax, theta)
    xi_tt = correlation_from_maps(temperature, temperature, lmax, theta)
    xi_ss = correlation_from_maps(smoothed, smoothed, lmax, theta)
    xi_ts = correlation_from_maps(temperature, smoothed, lmax, theta)
    xi_st = correlation_from_maps(smoothed, temperature, lmax, theta)
    trispectrum_projection = c_quadratic - xi_tt * xi_ss - xi_ts * xi_st
    return (
        correlation_from_maps(nonlinear_field, nonlinear_field, lmax, theta),
        correlation_from_maps(temperature, linear_field, lmax, theta),
        c_quadratic,
        trispectrum_projection,
    )


def sign_paired_descriptors(
    temperature: np.ndarray, lmax: int, smoothing_radians: float, linear_mark: float, quadratic_mark: float,
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Average t and -t to remove all odd-order contributions exactly."""
    plus = descriptors(temperature, lmax, smoothing_radians, linear_mark, quadratic_mark, theta)
    minus = descriptors(-temperature, lmax, smoothing_radians, linear_mark, quadratic_mark, theta)
    return tuple(0.5 * (positive + negative) for positive, negative in zip(plus, minus))


def correlation_to_cl(correlation: np.ndarray, mu: np.ndarray, weights: np.ndarray, lmax: int) -> np.ndarray:
    """Legendre-transform an angular correlation sampled at Gauss-Legendre nodes."""
    cl = np.empty(lmax + 1)
    previous = np.ones_like(mu)
    cl[0] = 2.0 * np.pi * np.sum(weights * correlation * previous)
    current = mu.copy()
    cl[1] = 2.0 * np.pi * np.sum(weights * correlation * current)
    for ell in range(2, lmax + 1):
        following = ((2.0 * ell - 1.0) * mu * current - (ell - 1.0) * previous) / ell
        cl[ell] = 2.0 * np.pi * np.sum(weights * correlation * following)
        previous, current = current, following
    return cl


def harmonic_descriptors(
    temperature: np.ndarray, lmax: int, smoothing_radians: float, linear_mark: float, quadratic_mark: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return C_ell descriptors, including the Wick-subtracted T_Q projection."""
    import healpy as hp

    temperature = temperature - temperature.mean()
    smoothed = smooth_map(temperature, lmax, smoothing_radians)
    nonlinear_field = temperature * (
        1.0 + linear_mark * smoothed + quadratic_mark * (smoothed**2 - np.mean(smoothed**2))
    )
    nonlinear_field -= nonlinear_field.mean()
    quadratic_field = temperature * smoothed
    quadratic_field -= quadratic_field.mean()
    nonlinear_cl = hp.anafast(nonlinear_field, lmax=lmax, iter=3)

    # Perform the Wick subtraction in configuration space, then transform the
    # connected angular correlation to C_ell.
    mu, weights = np.polynomial.legendre.leggauss(2 * lmax + 1)
    theta_nodes = np.arccos(mu)
    c_quadratic = correlation_from_maps(quadratic_field, quadratic_field, lmax, theta_nodes)
    xi_tt = correlation_from_maps(temperature, temperature, lmax, theta_nodes)
    xi_ss = correlation_from_maps(smoothed, smoothed, lmax, theta_nodes)
    xi_ts = correlation_from_maps(temperature, smoothed, lmax, theta_nodes)
    xi_st = correlation_from_maps(smoothed, temperature, lmax, theta_nodes)
    trispectrum_correlation = c_quadratic - xi_tt * xi_ss - xi_ts * xi_st
    return (
        nonlinear_cl,
        correlation_to_cl(c_quadratic, mu, weights, lmax),
        correlation_to_cl(trispectrum_correlation, mu, weights, lmax),
    )


def sign_paired_harmonic_descriptors(
    temperature: np.ndarray, lmax: int, smoothing_radians: float, linear_mark: float, quadratic_mark: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average harmonic descriptors of t and -t to remove odd-order terms."""
    plus = harmonic_descriptors(temperature, lmax, smoothing_radians, linear_mark, quadratic_mark)
    minus = harmonic_descriptors(-temperature, lmax, smoothing_radians, linear_mark, quadratic_mark)
    return tuple(0.5 * (positive + negative) for positive, negative in zip(plus, minus))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-map", type=Path, default=Path("archive/cmb_pilots/run_cmb_lensing_64/lensed_temperature.npy"))
    parser.add_argument("--lmax", type=int, default=128)
    parser.add_argument("--smoothing-degrees", type=float, default=2.0)
    parser.add_argument("--linear-mark", type=float, default=0.5)
    parser.add_argument("--quadratic-mark", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument("--output-dir", type=Path, default=Path("run_cmb_marked_64"))
    args = parser.parse_args()
    if args.lmax < 2 or args.smoothing_degrees <= 0:
        parser.error("lmax must be at least 2 and smoothing-degrees must be positive.")

    import healpy as hp

    temperature = np.load(args.input_map)
    nside = hp.get_nside(temperature)
    if args.lmax > 3 * nside - 1:
        parser.error("lmax cannot exceed 3*nside-1 for this map.")
    temperature = (temperature - temperature.mean()) / temperature.std()
    alm = hp.map2alm(temperature, lmax=args.lmax, iter=3)
    control_alm = phase_randomize_alm(alm, args.lmax, args.seed)
    control = hp.alm2map(control_alm, nside, lmax=args.lmax)
    control -= control.mean()

    theta = np.geomspace(np.deg2rad(0.25), np.pi, 180)
    smoothing_radians = np.deg2rad(args.smoothing_degrees)
    nonlinear_lensed, ebina_lensed, quadratic_lensed, trispectrum_lensed = sign_paired_descriptors(
        temperature, args.lmax, smoothing_radians, args.linear_mark, args.quadratic_mark, theta
    )
    nonlinear_control, ebina_control, quadratic_control, trispectrum_control = sign_paired_descriptors(
        control, args.lmax, smoothing_radians, args.linear_mark, args.quadratic_mark, theta
    )
    nonlinear_residual = nonlinear_lensed - nonlinear_control
    quadratic_residual = quadratic_lensed - quadratic_control
    trispectrum_residual = trispectrum_lensed - trispectrum_control
    wick_identity_residual = quadratic_residual - trispectrum_residual

    nonlinear_ell_lensed, quadratic_ell_lensed, trispectrum_ell_lensed = sign_paired_harmonic_descriptors(
        temperature, args.lmax, smoothing_radians, args.linear_mark, args.quadratic_mark
    )
    nonlinear_ell_control, quadratic_ell_control, trispectrum_ell_control = sign_paired_harmonic_descriptors(
        control, args.lmax, smoothing_radians, args.linear_mark, args.quadratic_mark
    )
    nonlinear_ell_residual = nonlinear_ell_lensed - nonlinear_ell_control
    quadratic_ell_residual = quadratic_ell_lensed - quadratic_ell_control
    trispectrum_ell_residual = trispectrum_ell_lensed - trispectrum_ell_control
    isolated_ell_component = args.linear_mark**2 * trispectrum_ell_residual
    remaining_ell_residual = nonlinear_ell_residual - isolated_ell_component

    ell = np.arange(args.lmax + 1)
    cl_lensed = hp.alm2cl(alm)
    cl_control = hp.alm2cl(control_alm)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.output_dir / "cmb_marked_statistics.csv",
        np.column_stack((np.rad2deg(theta), nonlinear_lensed, nonlinear_control, nonlinear_residual,
                         ebina_lensed, ebina_control, quadratic_lensed, quadratic_control, quadratic_residual,
                         trispectrum_lensed, trispectrum_control, trispectrum_residual, wick_identity_residual)),
        delimiter=",",
        header=("theta_degrees,C_nonlinear_lensed,C_nonlinear_phase_randomized,Delta_C_nonlinear,"
                "X_Ebina_lensed,X_Ebina_phase_randomized,C_Q_lensed,C_Q_phase_randomized,Delta_C_Q,"
                "T_Q_lensed,T_Q_phase_randomized,Delta_T_Q,Delta_C_Q_minus_Delta_T_Q"),
        comments="",
    )
    np.savetxt(args.output_dir / "power_spectra.csv", np.column_stack((ell, cl_lensed, cl_control)),
               delimiter=",", header="ell,C_ell_lensed,C_ell_phase_randomized", comments="")
    np.savetxt(
        args.output_dir / "cmb_trispectrum_decomposition_ell.csv",
        np.column_stack((ell, nonlinear_ell_lensed, nonlinear_ell_control, nonlinear_ell_residual,
                         quadratic_ell_lensed, quadratic_ell_control, quadratic_ell_residual,
                         trispectrum_ell_lensed, trispectrum_ell_control, trispectrum_ell_residual,
                         isolated_ell_component, remaining_ell_residual)),
        delimiter=",",
        header=("ell,C_ell_nonlinear_lensed,C_ell_nonlinear_phase_randomized,Delta_C_ell_nonlinear,"
                "C_ell_Q_lensed,C_ell_Q_phase_randomized,Delta_C_ell_Q,"
                "C_ell_T_Q_lensed,C_ell_T_Q_phase_randomized,Delta_C_ell_T_Q,"
                "a_squared_Delta_C_ell_T_Q,Delta_C_ell_nonlinear_minus_a_squared_Delta_C_ell_T_Q"),
        comments="",
    )

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Wrote CSV outputs; Matplotlib is unavailable for the diagnostic figure.")
    else:
        figure, axes = plt.subplots(1, 4, figsize=(20, 4.2))
        valid = ell >= 2
        axes[0].loglog(ell[valid], cl_lensed[valid], label="lensed CMB")
        axes[0].loglog(ell[valid], cl_control[valid], "--", label="phase-randomized")
        axes[0].set(xlabel=r"$\ell$", ylabel=r"$C_\ell$", title="Matched angular power spectrum")
        axes[0].legend()
        theta_degrees = np.rad2deg(theta)
        axes[1].semilogx(theta_degrees, nonlinear_lensed, color="C0", label="nonlinear marked (lensed)")
        axes[1].semilogx(theta_degrees, nonlinear_control, "--", color="C0", label="nonlinear marked (control)")
        axes[1].semilogx(theta_degrees, ebina_lensed, color="C2", label="linear mark / Ebina (lensed)")
        axes[1].semilogx(theta_degrees, ebina_control, "--", color="C2", label="linear mark / Ebina (control)")
        axes[1].set(xlabel=r"$\theta$ [deg]", ylabel="angular correlation", title="CMB marked statistics")
        axes[1].legend(fontsize="small")
        axes[2].loglog(theta_degrees, np.abs(nonlinear_residual), color="C3", label="nonlinear marked residual")
        axes[2].set(xlabel=r"$\theta$ [deg]", ylabel="residual magnitude", title="Lensing-trispectrum response")
        axes[2].legend()
        axes[3].loglog(theta_degrees, np.abs(quadratic_residual), color="C4", label=r"$|\Delta C_Q|$")
        axes[3].loglog(theta_degrees, np.abs(trispectrum_residual), "--", color="C1", label=r"$|\Delta T_Q|$")
        axes[3].loglog(theta_degrees, np.abs(wick_identity_residual), ":", color="0.35",
                       label=r"$|\Delta C_Q-\Delta T_Q|$")
        axes[3].set(xlabel=r"$\theta$ [deg]", ylabel="magnitude", title="Direct trispectrum projection")
        axes[3].legend(fontsize="small")
        figure.tight_layout()
        figure.savefig(args.output_dir / "cmb_marked_statistics.png", dpi=180)

        valid = ell >= 2
        decomposition, decomposition_axes = plt.subplots(1, 2, figsize=(10, 4.2))
        decomposition_axes[0].axhline(0.0, color="0.55", lw=1)
        decomposition_axes[0].semilogx(ell[valid], nonlinear_ell_residual[valid], color="C3", label=r"$\Delta C_\ell^{\rm NL}$")
        decomposition_axes[0].semilogx(
            ell[valid], isolated_ell_component[valid], color="C1", ls="--", label=r"$a^2\Delta C_{\ell,T_Q}$"
        )
        decomposition_axes[0].set(xlabel=r"$\ell$", ylabel="signed residual", title="Isolated trispectrum contribution")
        decomposition_axes[0].legend(fontsize="small")
        decomposition_axes[1].loglog(ell[valid], np.abs(nonlinear_ell_residual[valid]), color="C3", label=r"$|\Delta C_\ell^{\rm NL}|$")
        decomposition_axes[1].loglog(
            ell[valid], np.abs(isolated_ell_component[valid]), color="C1", ls="--", label=r"$|a^2\Delta C_{\ell,T_Q}|$"
        )
        decomposition_axes[1].loglog(
            ell[valid], np.abs(remaining_ell_residual[valid]), color="C4", ls=":", label="remaining response"
        )
        decomposition_axes[1].set(xlabel=r"$\ell$", ylabel="magnitude", title=r"$\Delta C_\ell^{\rm NL}-a^2\Delta C_{\ell,T_Q}$")
        decomposition_axes[1].legend(fontsize="small")
        decomposition.tight_layout()
        decomposition.savefig(args.output_dir / "cmb_trispectrum_decomposition_ell.png", dpi=180)

    max_relative_cl_difference = np.max(np.abs(cl_lensed[2:] - cl_control[2:]) / np.maximum(cl_lensed[2:], 1e-30))
    print(f"Wrote results to {args.output_dir}")
    print(f"nside / lmax: {nside} / {args.lmax}")
    print(f"max relative C_ell difference: {max_relative_cl_difference:.3e}")
    print(f"max |nonlinear marked residual|: {np.max(np.abs(nonlinear_residual)):.3e}")
    print(f"max |Delta C_Q - Delta T_Q|: {np.max(np.abs(wick_identity_residual)):.3e}")


if __name__ == "__main__":
    main()
