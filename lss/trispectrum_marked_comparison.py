#!/usr/bin/env python3
"""Compare Ebina's linear mark with a nonlinear marked auto-correlation.

This is a B=0, T!=0 toy experiment.  The non-Gaussian field is an odd local
cubic transform of a Gaussian field.  Results are paired between delta and
-delta, so every odd (including bispectrum) contribution cancels exactly,
whereas the connected four-point contribution remains.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from lss.lognormal_marked_correlation import (
    correlation_grid,
    gaussian_random_field,
    isotropic_power,
    phase_randomize,
    radial_average,
    smooth_field,
)


def cubic_symmetric_field(gaussian: np.ndarray, amplitude: float) -> np.ndarray:
    """Return an odd local transform with vanishing ensemble bispectrum.

    ``g**3 - 3 <g**2> g`` is the third Hermite operator.  It changes the
    connected four-point function at leading order while preserving the global
    sign symmetry that forces all odd moments to vanish.
    """
    variance = gaussian.var()
    field = gaussian + amplitude * (gaussian**3 - 3.0 * variance * gaussian)
    return field - field.mean()


def cross_correlation_grid(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Periodic cross-correlation <left(x) right(x+r)> on the mesh."""
    left_k = np.fft.rfftn(left)
    right_k = np.fft.rfftn(right)
    return np.fft.irfftn(left_k.conj() * right_k, s=left.shape).real / left.size


def normalized_mark(delta: np.ndarray, boxsize: float, radius: float, power: float) -> np.ndarray:
    """Return F/<F>-1 for F=(1+delta)(1+delta_R)^power."""
    smoothed = smooth_field(delta, boxsize, radius)
    weighted = (1.0 + delta) * (1.0 + smoothed) ** power
    return weighted / weighted.mean() - 1.0


def ebina_linear_mark(delta: np.ndarray, boxsize: float, radius: float) -> np.ndarray:
    """Return the density-weighted linear marked fluctuation, C0=C1=1."""
    smoothed = smooth_field(delta, boxsize, radius)
    weighted = (1.0 + delta) * (1.0 + smoothed)
    return weighted / weighted.mean() - 1.0


def quadratic_trispectrum_projection(
    delta: np.ndarray, boxsize: float, radius: float, nbins: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Measure the collapsed connected four-point projection of ``delta*delta_R``.

    For Q=delta*s-<delta*s>, C_Q=<Q_1 Q_2> contains one connected
    trispectrum term plus the two non-local Wick contractions.  Subtracting
    those contractions gives T_Q exactly.
    """
    smoothed = smooth_field(delta, boxsize, radius)
    quadratic = delta * smoothed
    quadratic -= quadratic.mean()
    radii, c_quadratic = radial_average(correlation_grid(quadratic), boxsize, nbins)
    _, xi_dd = radial_average(correlation_grid(delta), boxsize, nbins)
    _, xi_ss = radial_average(correlation_grid(smoothed), boxsize, nbins)
    _, xi_ds = radial_average(cross_correlation_grid(delta, smoothed), boxsize, nbins)
    _, xi_sd = radial_average(cross_correlation_grid(smoothed, delta), boxsize, nbins)
    connected = c_quadratic - xi_dd * xi_ss - xi_ds * xi_sd
    return radii, c_quadratic, connected


def paired_statistics(
    delta: np.ndarray, control: np.ndarray, boxsize: float, radius: float, mark_power: float, nbins: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Average statistics over sign partners for exact odd-order cancellation."""
    marked_ng = []
    marked_control = []
    ebina_ng = []
    ebina_control = []
    radii = None
    for field, destination, cross_destination in (
        (delta, marked_ng, ebina_ng),
        (-delta, marked_ng, ebina_ng),
        (control, marked_control, ebina_control),
        (-control, marked_control, ebina_control),
    ):
        marked = normalized_mark(field, boxsize, radius, mark_power)
        ebina = ebina_linear_mark(field, boxsize, radius)
        radii, auto = radial_average(correlation_grid(marked), boxsize, nbins)
        _, cross = radial_average(cross_correlation_grid(field, ebina), boxsize, nbins)
        destination.append(1.0 + auto)
        cross_destination.append(cross)
    return (
        radii,
        np.mean(marked_ng, axis=0),
        np.mean(marked_control, axis=0),
        np.mean(ebina_ng, axis=0),
        np.mean(ebina_control, axis=0),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmesh", type=int, default=64)
    parser.add_argument("--boxsize", type=float, default=500.0)
    parser.add_argument("--sigma-g", type=float, default=0.25)
    parser.add_argument("--cubic-amplitude", type=float, default=2.0)
    parser.add_argument("--k0", type=float, default=0.08)
    parser.add_argument("--slope", type=float, default=2.0)
    parser.add_argument("--smoothing-radius", type=float, default=10.0)
    parser.add_argument("--mark-power", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--nbins", type=int, default=50)
    parser.add_argument("--output-dir", type=Path, default=Path("run_trispectrum_64"))
    args = parser.parse_args()
    if args.nmesh < 4 or args.nmesh % 2:
        parser.error("--nmesh must be an even integer of at least 4.")
    if args.boxsize <= 0 or args.smoothing_radius <= 0 or args.sigma_g <= 0:
        parser.error("box size, smoothing radius, and sigma-g must be positive.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    gaussian = gaussian_random_field(args.nmesh, args.boxsize, args.sigma_g, args.seed, args.k0, args.slope)
    delta = cubic_symmetric_field(gaussian, args.cubic_amplitude)
    control = phase_randomize(delta, args.seed + 1)
    control -= control.mean() - delta.mean()

    radii, marked_ng, marked_control, ebina_ng, ebina_control = paired_statistics(
        delta, control, args.boxsize, args.smoothing_radius, args.mark_power, args.nbins
    )
    auto_residual = marked_ng - marked_control
    ebina_residual = ebina_ng - ebina_control
    _, quadratic_ng, trispectrum_ng = quadratic_trispectrum_projection(
        delta, args.boxsize, args.smoothing_radius, args.nbins
    )
    _, quadratic_control, trispectrum_control = quadratic_trispectrum_projection(
        control, args.boxsize, args.smoothing_radius, args.nbins
    )
    quadratic_residual = quadratic_ng - quadratic_control
    trispectrum_residual = trispectrum_ng - trispectrum_control
    wick_identity_residual = quadratic_residual - trispectrum_residual
    isolated_mark_component = args.mark_power**2 * trispectrum_residual
    remaining_mark_residual = auto_residual - isolated_mark_component
    k, power_ng = isotropic_power(delta, args.boxsize, args.nbins)
    _, power_control = isotropic_power(control, args.boxsize, args.nbins)

    np.savetxt(
        args.output_dir / "trispectrum_comparison.csv",
        np.column_stack((radii, marked_ng, marked_control, auto_residual, ebina_ng, ebina_control, ebina_residual,
                         quadratic_ng, quadratic_control, quadratic_residual, trispectrum_ng, trispectrum_control,
                         trispectrum_residual, wick_identity_residual, isolated_mark_component, remaining_mark_residual)),
        delimiter=",",
        header=("r,M_nonlinear_NG,M_nonlinear_phase_randomized,Delta_M,"
                "X_Ebina_NG,X_Ebina_phase_randomized,Delta_X_Ebina,"
                "C_Q_NG,C_Q_phase_randomized,Delta_C_Q,T_Q_NG,T_Q_phase_randomized,Delta_T_Q,"
                "Delta_C_Q_minus_Delta_T_Q,p_squared_Delta_T_Q,Delta_M_minus_p_squared_Delta_T_Q"),
        comments="",
    )
    np.savetxt(
        args.output_dir / "power_spectra.csv", np.column_stack((k, power_ng, power_control)),
        delimiter=",", header="k,P_cubic_symmetric,P_phase_randomized", comments=""
    )

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed; wrote CSV outputs but skipped the diagnostic plot.")
    else:
        figure, axes = plt.subplots(1, 4, figsize=(20, 4.2))
        axes[0].plot(k, power_ng, label=r"$B=0,\ T\ne0$")
        axes[0].plot(k, power_control, "--", label="phase-randomized")
        axes[0].set(xlabel=r"$k$", ylabel=r"$P(k)$", yscale="log", title="Matched power spectrum")
        axes[0].legend()
        axes[1].plot(radii, marked_ng - 1.0, color="C0", label=r"nonlinear $M-1$ ($B=0,\ T\ne0$)")
        axes[1].plot(radii, marked_control - 1.0, "--", color="C0", label=r"nonlinear $M-1$ (phase-randomized)")
        axes[1].plot(radii, ebina_ng, color="C2", label=r"Ebina $X$ ($B=0,\ T\ne0$)")
        axes[1].plot(radii, ebina_control, "--", color="C2", label=r"Ebina $X$ (phase-randomized)")
        axes[1].set(xlabel=r"$r$", ylabel="correlation excess", xscale="log",
                    title="Nonlinear and Ebina marked statistics")
        axes[1].legend(fontsize="small")
        axes[2].loglog(radii, np.abs(auto_residual), color="C3", label=r"$|\Delta M|$ (nonlinear auto)")
        axes[2].set(xlabel=r"$r$", ylabel="residual magnitude", title="Trispectrum-only response")
        axes[2].legend(loc="upper right", fontsize="small")
        axes[3].loglog(radii, np.abs(quadratic_residual), color="C4", label=r"$|\Delta C_Q|$")
        axes[3].loglog(radii, np.abs(trispectrum_residual), "--", color="C1", label=r"$|\Delta T_Q|$")
        axes[3].loglog(radii, np.abs(wick_identity_residual), ":", color="0.35",
                       label=r"$|\Delta C_Q-\Delta T_Q|$")
        axes[3].set(xlabel=r"$r$", ylabel="magnitude", title="Direct trispectrum projection")
        axes[3].legend(fontsize="small")
        figure.tight_layout()
        figure.savefig(args.output_dir / "trispectrum_comparison.png", dpi=180)

        decomposition, decomposition_axes = plt.subplots(1, 2, figsize=(10, 4.2))
        decomposition_axes[0].axhline(0.0, color="0.55", lw=1)
        decomposition_axes[0].semilogx(radii, auto_residual, color="C3", label=r"$\Delta M$")
        decomposition_axes[0].semilogx(
            radii, isolated_mark_component, color="C1", ls="--", label=r"$p^2\Delta T_{\delta s,\delta s}$"
        )
        decomposition_axes[0].set(
            xlabel=r"$r$", ylabel="signed residual", title="Isolated trispectrum contribution"
        )
        decomposition_axes[0].legend(fontsize="small")
        decomposition_axes[1].loglog(radii, np.abs(auto_residual), color="C3", label=r"$|\Delta M|$")
        decomposition_axes[1].loglog(
            radii, np.abs(isolated_mark_component), color="C1", ls="--", label=r"$|p^2\Delta T_{\delta s,\delta s}|$"
        )
        decomposition_axes[1].loglog(
            radii, np.abs(remaining_mark_residual), color="C4", ls=":", label="remaining response"
        )
        decomposition_axes[1].set(xlabel=r"$r$", ylabel="magnitude", title=r"$\Delta M-p^2\Delta T_{\delta s,\delta s}$")
        decomposition_axes[1].legend(fontsize="small")
        decomposition.tight_layout()
        decomposition.savefig(args.output_dir / "trispectrum_decomposition.png", dpi=180)

    print(f"Wrote results to {args.output_dir}")
    print(f"max relative power difference: {np.nanmax(np.abs(power_ng - power_control) / np.maximum(power_ng, 1e-30)):.3e}")
    print(f"max |nonlinear auto residual|: {np.nanmax(np.abs(auto_residual)):.3e}")
    print(f"max |Ebina cross residual|: {np.nanmax(np.abs(ebina_residual)):.3e}")
    print(f"max |Delta C_Q - Delta T_Q|: {np.nanmax(np.abs(wick_identity_residual)):.3e}")
    print(f"max |remaining marked residual|: {np.nanmax(np.abs(remaining_mark_residual)):.3e}")


if __name__ == "__main__":
    main()
