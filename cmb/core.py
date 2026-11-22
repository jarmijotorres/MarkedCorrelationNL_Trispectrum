"""Shared CMB spectra, phase controls, and angular projection."""

from __future__ import annotations

import numpy as np


def bin_spectrum(cl: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Average C_L in half-open bins with 2L+1 mode weights."""
    ell = np.arange(cl.size)
    return np.array([
        np.average(cl[lo:hi], weights=2 * ell[lo:hi] + 1)
        for lo, hi in zip(edges[:-1], edges[1:])
    ])


def toy_spectra(lmax: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return spectra for a low-resolution software check."""
    ell = np.arange(lmax + 1)
    tt = np.zeros(lmax + 1)
    pp = np.zeros(lmax + 1)
    valid = ell >= 2
    tt[valid] = 2.0e-5 / (ell[valid] * (ell[valid] + 1.0))
    pp[valid] = 5.0e-7 / (ell[valid] * (ell[valid] + 1.0)) ** 2
    return ell, tt, pp


def phase_randomize_alm(alm: np.ndarray, lmax: int, seed: int) -> np.ndarray:
    """Randomize HEALPix alm phases while preserving every mode amplitude."""
    import healpy as hp

    rng = np.random.default_rng(seed)
    ell, emm = hp.Alm.getlm(lmax)
    phases = np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, size=alm.size))
    phases[emm == 0] = rng.choice(np.array([-1.0, 1.0]), size=np.count_nonzero(emm == 0))
    phases[ell == 0] = 1.0
    return alm * phases


def correlation_to_cl(
    correlation: np.ndarray, mu: np.ndarray, weights: np.ndarray, lmax: int,
) -> np.ndarray:
    """Project an angular correlation at Gauss-Legendre nodes to C_ell."""
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
