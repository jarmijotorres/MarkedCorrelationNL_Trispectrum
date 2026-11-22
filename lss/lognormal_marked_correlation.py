#!/usr/bin/env python3
"""A harmonic-space toy model for nonlinear marked correlations.

Generate a periodic lognormal density field and a phase-randomized control with
identical Fourier amplitudes.  The difference between their marked
correlations isolates the higher-order spatial information of the lognormal
realization while holding its two-point function fixed.

Example
-------
python3 -m lss.lognormal_marked_correlation --nmesh 128 --boxsize 500 \
    --smoothing-radius 10 --mark-power 2 --output-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def k_grids(nmesh: int, boxsize: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return Fourier coordinate arrays compatible with ``numpy.rfftn``."""
    spacing = boxsize / nmesh
    kx = 2.0 * np.pi * np.fft.fftfreq(nmesh, d=spacing)
    ky = 2.0 * np.pi * np.fft.fftfreq(nmesh, d=spacing)
    kz = 2.0 * np.pi * np.fft.rfftfreq(nmesh, d=spacing)
    kx3, ky3, kz3 = np.meshgrid(kx, ky, kz, indexing="ij")
    return kx3, ky3, kz3, np.sqrt(kx3**2 + ky3**2 + kz3**2)


def input_power_spectrum(k: np.ndarray, amplitude: float, k0: float, slope: float) -> np.ndarray:
    """Smooth toy spectrum for the Gaussian log-density field.

    It is finite in the infrared and exponentially suppressed in the UV;
    ``amplitude`` is subsequently adjusted through the measured field variance.
    """
    return amplitude * (1.0 + (k / k0) ** 2) ** (-0.5 * slope) * np.exp(-(k / (4.0 * k0)) ** 2)


def gaussian_random_field(
    nmesh: int, boxsize: float, sigma_g: float, seed: int, k0: float, slope: float
) -> np.ndarray:
    """Draw a real Gaussian field with approximately the requested variance."""
    rng = np.random.default_rng(seed)
    _, _, _, kmag = k_grids(nmesh, boxsize)
    white = rng.normal(size=(nmesh, nmesh, nmesh))
    white_k = np.fft.rfftn(white)
    shape = input_power_spectrum(kmag, 1.0, k0, slope)
    shape[0, 0, 0] = 0.0
    field = np.fft.irfftn(white_k * np.sqrt(shape), s=(nmesh,) * 3).real
    field -= field.mean()
    field *= sigma_g / field.std()
    return field


def phase_randomize(field: np.ndarray, seed: int) -> np.ndarray:
    """Return a real field with exactly the same discrete Fourier amplitudes.

    A real white-noise field supplies phases with the required Hermitian
    symmetry, including the rFFT boundary planes.  Consequently, the inverse
    transform is real and the input and output auto-power agree mode by mode.
    """
    rng = np.random.default_rng(seed)
    field_k = np.fft.rfftn(field)
    phase_source = np.fft.rfftn(rng.normal(size=field.shape))
    phase = np.ones_like(phase_source)
    nonzero = np.abs(phase_source) > 0.0
    phase[nonzero] = phase_source[nonzero] / np.abs(phase_source[nonzero])
    randomized_k = np.abs(field_k) * phase
    randomized_k[0, 0, 0] = field_k[0, 0, 0]
    return np.fft.irfftn(randomized_k, s=field.shape).real


def smooth_field(field: np.ndarray, boxsize: float, radius: float) -> np.ndarray:
    """Gaussian smooth a periodic field, W_R(k)=exp(-k^2 R^2/2)."""
    _, _, _, kmag = k_grids(field.shape[0], boxsize)
    window = np.exp(-0.5 * (kmag * radius) ** 2)
    return np.fft.irfftn(np.fft.rfftn(field) * window, s=field.shape).real


def marked_fluctuation(delta: np.ndarray, boxsize: float, radius: float, mark_power: float) -> np.ndarray:
    """Build F/<F>-1 for F=(1+delta)(1+delta_R)^p."""
    smoothed = smooth_field(delta, boxsize, radius)
    base = 1.0 + smoothed
    if not float(mark_power).is_integer() and np.any(base <= 0.0):
        raise ValueError("1 + delta_R must be positive for a non-integer mark power.")
    weighted = (1.0 + delta) * base**mark_power
    return weighted / weighted.mean() - 1.0


def correlation_grid(field: np.ndarray) -> np.ndarray:
    """Periodic auto-correlation <f(x)f(x+r)> evaluated on the mesh."""
    fourier = np.fft.rfftn(field)
    return np.fft.irfftn(np.abs(fourier) ** 2, s=field.shape).real / field.size


def radial_average(grid: np.ndarray, boxsize: float, nbins: int) -> tuple[np.ndarray, np.ndarray]:
    """Spherically average a periodic-grid statistic up to half the box size."""
    nmesh = grid.shape[0]
    spacing = boxsize / nmesh
    coordinate = np.minimum(np.arange(nmesh), nmesh - np.arange(nmesh)) * spacing
    x, y, z = np.meshgrid(coordinate, coordinate, coordinate, indexing="ij")
    radius = np.sqrt(x**2 + y**2 + z**2).ravel()
    values = grid.ravel()
    edges = np.linspace(0.0, boxsize / 2.0, nbins + 1)
    indices = np.digitize(radius, edges) - 1
    valid = (indices >= 0) & (indices < nbins)
    count = np.bincount(indices[valid], minlength=nbins)
    total = np.bincount(indices[valid], weights=values[valid], minlength=nbins)
    average = np.divide(total, count, out=np.full(nbins, np.nan), where=count > 0)
    return 0.5 * (edges[1:] + edges[:-1]), average


def isotropic_power(field: np.ndarray, boxsize: float, nbins: int) -> tuple[np.ndarray, np.ndarray]:
    """Return a shell-averaged power estimate; normalization is P in box units."""
    nmesh = field.shape[0]
    _, _, _, kmag = k_grids(nmesh, boxsize)
    modes = np.fft.rfftn(field)
    power = boxsize**3 * np.abs(modes) ** 2 / nmesh**6
    nonzero_k = kmag[kmag > 0]
    edges = np.linspace(nonzero_k.min(), nonzero_k.max(), nbins + 1)
    indices = np.digitize(kmag.ravel(), edges) - 1
    valid = (indices >= 0) & (indices < nbins)
    count = np.bincount(indices[valid], minlength=nbins)
    total = np.bincount(indices[valid], weights=power.ravel()[valid], minlength=nbins)
    average = np.divide(total, count, out=np.full(nbins, np.nan), where=count > 0)
    return 0.5 * (edges[1:] + edges[:-1]), average


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nmesh", type=int, default=128, help="Cubic mesh size (default: 128).")
    parser.add_argument("--boxsize", type=float, default=500.0, help="Periodic box side length.")
    parser.add_argument("--sigma-g", type=float, default=0.7, help="Standard deviation of the Gaussian log-density.")
    parser.add_argument("--k0", type=float, default=0.08, help="Turnover scale of the input Gaussian spectrum.")
    parser.add_argument("--slope", type=float, default=2.0, help="High-k slope of the input Gaussian spectrum.")
    parser.add_argument("--smoothing-radius", type=float, default=10.0, help="Gaussian smoothing radius R.")
    parser.add_argument("--mark-power", type=float, default=2.0, help="Exponent p of (1+delta_R)^p.")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for the lognormal realization.")
    parser.add_argument("--nbins", type=int, default=50, help="Number of radial bins.")
    parser.add_argument("--output-dir", type=Path, default=Path("lognormal_marked_output"))
    args = parser.parse_args()
    if args.nmesh < 4 or args.nmesh % 2:
        parser.error("--nmesh must be an even integer of at least 4.")
    if args.boxsize <= 0 or args.smoothing_radius <= 0 or args.sigma_g <= 0:
        parser.error("box size, smoothing radius, and sigma-g must be positive.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    gaussian = gaussian_random_field(args.nmesh, args.boxsize, args.sigma_g, args.seed, args.k0, args.slope)
    delta_ln = np.exp(gaussian - 0.5 * gaussian.var()) - 1.0
    delta_control = phase_randomize(delta_ln, args.seed + 1)
    delta_control -= delta_control.mean() - delta_ln.mean()

    marked_ln = marked_fluctuation(delta_ln, args.boxsize, args.smoothing_radius, args.mark_power)
    marked_control = marked_fluctuation(delta_control, args.boxsize, args.smoothing_radius, args.mark_power)
    radii, xi_ln = radial_average(correlation_grid(delta_ln), args.boxsize, args.nbins)
    _, xi_control = radial_average(correlation_grid(delta_control), args.boxsize, args.nbins)
    _, xi_mark_ln = radial_average(correlation_grid(marked_ln), args.boxsize, args.nbins)
    _, xi_mark_control = radial_average(correlation_grid(marked_control), args.boxsize, args.nbins)
    marked_ratio_ln = (1.0 + xi_mark_ln) / (1.0 + xi_ln)
    marked_ratio_control = (1.0 + xi_mark_control) / (1.0 + xi_control)
    delta_marked = xi_mark_ln - xi_mark_control

    k, power_ln = isotropic_power(delta_ln, args.boxsize, args.nbins)
    _, power_control = isotropic_power(delta_control, args.boxsize, args.nbins)
    np.savetxt(
        args.output_dir / "marked_correlations.csv",
        np.column_stack((radii, xi_ln, xi_control, 1 + xi_mark_ln, 1 + xi_mark_control, marked_ratio_ln, marked_ratio_control, delta_marked)),
        delimiter=",",
        header="r,xi_lognormal,xi_phase_randomized,M_lognormal,M_phase_randomized,marked_ratio_lognormal,marked_ratio_phase_randomized,Delta_M",
        comments="",
    )
    np.savetxt(args.output_dir / "power_spectra.csv", np.column_stack((k, power_ln, power_control)), delimiter=",", header="k,P_lognormal,P_phase_randomized", comments="")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed; wrote CSV outputs but skipped the diagnostic plot.")
    else:
        figure, axes = plt.subplots(1, 3, figsize=(15, 4.2))
        axes[0].plot(k, power_ln, label="lognormal")
        axes[0].plot(k, power_control, "--", label="phase-randomized")
        axes[0].set(xlabel=r"$k$", ylabel=r"$P(k)$", yscale="log", title="Matched power spectrum")
        axes[0].legend()
        axes[1].plot(radii, 1 + xi_mark_ln, label="lognormal")
        axes[1].plot(radii, 1 + xi_mark_control, "--", label="phase-randomized")
        axes[1].set(xlabel=r"$r$", ylabel=r"$M(r)$", title="Marked auto-correlation")
        axes[1].legend()
        axes[2].loglog(radii, np.abs(delta_marked), color="C3")
        axes[2].set(
            xlabel=r"$r$",
            ylabel=r"$|\Delta M(r)|$",
            title="Higher-order residual magnitude",
        )
        figure.tight_layout()
        figure.savefig(args.output_dir / "marked_correlation_diagnostics.png", dpi=180)

    max_relative_power_difference = np.nanmax(np.abs(power_ln - power_control) / np.maximum(power_ln, 1e-30))
    print(f"Wrote results to {args.output_dir}")
    print(f"delta mean / variance: {delta_ln.mean():.3e} / {delta_ln.var():.3e}")
    print(f"maximum binned relative power difference: {max_relative_power_difference:.3e}")
    print(f"max |Delta M| (nonzero r bins): {np.nanmax(np.abs(delta_marked[1:])):.3e}")


if __name__ == "__main__":
    main()
