#!/usr/bin/env python3
"""Generate a simple full-sky CMB lensing trispectrum realization.

The unlensed CMB temperature and lensing potential are independent Gaussian
fields (C_ell^{T phi}=0).  Lensing then generates a connected trispectrum,
while the ensemble temperature bispectrum remains zero.  The spectra below
are deliberately simple toy spectra for an integration test, not a CAMB
cosmology.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def toy_spectra(lmax: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ell, C_ell^TT and C_ell^phiphi for a low-resolution smoke test."""
    ell = np.arange(lmax + 1)
    tt = np.zeros(lmax + 1)
    pp = np.zeros(lmax + 1)
    valid = ell >= 2
    tt[valid] = 2.0e-5 / (ell[valid] * (ell[valid] + 1.0))
    pp[valid] = 5.0e-7 / (ell[valid] * (ell[valid] + 1.0)) ** 2
    return ell, tt, pp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--lmax", type=int, default=128)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--output-dir", type=Path, default=Path("run_cmb_lensing_64"))
    args = parser.parse_args()
    if args.nside < 1 or args.lmax < 2:
        parser.error("nside must be positive and lmax must be at least 2.")

    from lenspyx.lensing import synfast

    ell, tt, pp = toy_spectra(args.lmax)
    maps = synfast(
        {"tt": tt, "pp": pp},  # no ``tp`` term: primary T and phi are independent
        lmax=args.lmax,
        geometry=("healpix", {"nside": args.nside}),
        seed=args.seed,
        verbose=False,
    )
    temperature = np.asarray(maps["T"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.save(args.output_dir / "lensed_temperature.npy", temperature)
    np.savetxt(
        args.output_dir / "input_spectra.csv",
        np.column_stack((ell, tt, pp)),
        delimiter=",",
        header="ell,C_ell_TT,C_ell_phiphi",
        comments="",
    )

    try:
        import healpy as hp
        import matplotlib.pyplot as plt
    except ImportError:
        print("Wrote the map and spectra; install healpy and matplotlib to make a Mollweide image.")
    else:
        hp.mollview(temperature, title="Toy lensed CMB temperature ($T$ independent of $\\phi$)", unit="arb. units")
        plt.savefig(args.output_dir / "lensed_temperature_mollweide.png", dpi=180, bbox_inches="tight")
        plt.close()

    print(f"Wrote results to {args.output_dir}")
    print(f"pixels: {temperature.size}; mean: {temperature.mean():.3e}; standard deviation: {temperature.std():.3e}")


if __name__ == "__main__":
    main()
