# Current CMB marked-correlation configuration

## Simulation and temperature maps

- Generate an independent Gaussian primary CMB temperature field \(T\) and Gaussian lensing potential \(\phi\) from fiducial CAMB spectra. Lens the primary with \(\phi\) using lenspyx. Keep the same primary as the paired unlensed map.
- Use full-sky HEALPix maps with \(N_{\rm side}=2048\). Generate input harmonics through \(\ell_{\max}=3000\).
- For **each** lensed or unlensed temperature map, retain \(100\leq\ell\leq3000\), subtract its mean, and divide by its own standard deviation. Call the resulting dimensionless contrast \(t\).

## Smoothing and mark

Smooth the filtered, normalized temperature contrast separately in each case:

\[
s_{{\rm len},\ell m}=W_\ell t_{{\rm len},\ell m},\qquad
s_{{\rm unlen},\ell m}=W_\ell t_{{\rm unlen},\ell m},\qquad
W_\ell=\exp\!\left[-\tfrac12\ell(\ell+1)\sigma_s^2\right].
\]

The Gaussian width is **\(\sigma_s=2\) arcmin**, or FWHM \(4.71\) arcmin. The same window is used for both temperature maps. The potential \(\phi\) supplies the lensing remapping and the \(q\)–\(\phi\) diagnostic. It is not smoothed to construct \(s\).

Use the fixed \(p=2\) mark and its two sign partners:

\[
F_\pm=(1\pm t)(1\pm s)^2,\qquad
\mu_\pm=\langle F_\pm\rangle,\qquad
\delta_{m,\pm}=F_\pm/\mu_\pm-1.
\]

Each sign partner uses its own mean \(\mu_\pm\). Sign pairing averages their auto-spectra to remove odd-order terms in \(t\).

## Measured statistics

The primary marked spectrum is

\[
C_L^{\rm mark}=\tfrac12\left(C_L^{\delta_{m,+}\delta_{m,+}}+
C_L^{\delta_{m,-}\delta_{m,-}}\right).
\]

Its angular counterpart is \(M(\theta)=1+\xi_{\delta_m\delta_m}(\theta)\), with the same sign averaging understood. Here \(\ell\) labels input temperature modes and \(L\) labels output modes of the marked and quadratic fields.

For the explicit connected four-point projection, form \(q=ts-\langle ts\rangle\) and calculate

\[
T_{q,L}=C_L^{qq}-W_L^{qq},\qquad
W^{qq}(\theta)=\xi_{tt}(\theta)\xi_{ss}(\theta)+\xi_{ts}(\theta)^2.
\]

The identified \(q\) contribution in the sign-averaged marked spectrum has coefficient \(K_q=2(\mu_+^{-2}+\mu_-^{-2})\). Compare \(K_qT_{q,L}\) with the full marked spectrum and its remainder. Also measure \(C_L^{q\phi}\) using the known input potential as a lensing-response diagnostic.

For each temperature map, construct phase-randomized controls that preserve its measured \(C_\ell^{tt}\). Subtract the mean control statistic for each lensed and unlensed map, then take their paired difference. The current pilot has **one realization and one control per map**. The intended ensemble test will use multiple independent paired seeds and multiple controls per map.

## Output range and resolution check

The current calculation saves unbinned spectra over \(40\leq L<3000\). Its plotting bins have width \(\Delta L=10\) below \(L=200\), 25 from \(L=200\) to 999, and 50 from \(L=1000\) to 2999. Short-angle results should be formed from the saved unbinned spectra, retaining the exact Legendre weights of each \(L\).

The connected response observed in this pilot is concentrated mainly around output \(100\leq L<600\) at separations \(2'\)–\(10'\). The present \(2'\) angular value uses output modes only through \(L=2999\). Products of input temperature modes through \(\ell=3000\) can populate \(q\) up to roughly \(L=6000\) and the cubic term in the marked field up to roughly \(L=9000\). The higher output modes and pixel resolution require a convergence check before interpreting the full small-angle marked correlation. The low-\(L\) response and the small-angle correlation should both be compared on matched seeds at a finer resolution.

The \(\ell\geq100\) temperature cut and \(\sigma_s=2'\) smoothing define the saved angular-consistency run. The \(\ell\geq600\) run provides a matched filter comparison. That single-realization comparison did not establish an optimal filter, smoothing scale, or detection significance.

## Subsequent marked-spectrum screen

A later matched three-seed pilot with three phase controls per map compared the \(\ell\geq100\) and \(\ell\geq600\) input cuts at fixed smoothing and mark. The paired marked \(C_L\), averaged over output \(100\leq L<600\), changed sign across the three \(\ell\geq100\) seeds. It was positive in all three \(\ell\geq600\) seeds, with bandpowers \(6.19\), \(7.09\), and \(6.97\times10^{-8}\). The full marked spectrum is the target statistic. The \(q\) component is retained as a check of its four-point contribution. See [the pilot analysis](run_cmb_marked_pilot_lmin600_2048/analysis.md). The \(\ell\geq600\) input is now the leading candidate for the next ensemble, pending control-variance and resolution checks.
