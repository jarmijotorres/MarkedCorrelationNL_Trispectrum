# CMB lensing marked correlation pilot

## Current working configuration

The current baseline uses one full-sky, noiseless, independent Gaussian primary temperature and Gaussian lensing-potential realization drawn from CAMB spectra. Lenspyx lenses the primary with the generated potential. The same primary supplies the paired unlensed map. The map resolution is `Nside=2048`, and the input harmonic calculation ends at `ell=3000`.

For each lensed or unlensed map, retain input temperature modes `100 <= ell <= 3000`. Center that filtered map and divide by its own standard deviation to obtain the dimensionless temperature contrast `t`. The resulting amplitudes are specific to this variance-normalized construction. The earlier `ell >= 600` run is a comparison with the same realization, not a different smoothing scale.

Smooth **the filtered, normalized temperature contrast** with a Gaussian harmonic window,

\[
s_{\ell m}=W_\ell t_{\ell m},\qquad
W_\ell=\exp[-\ell(\ell+1)\sigma_s^2/2],\qquad
\sigma_s=2\,\mathrm{arcmin}.
\]

The smoothing FWHM is 4.71 arcmin. The window is 0.998 at `ell=100`, 0.941 at 600, 0.508 at 2000, and 0.218 at 3000. The lensing potential is used for remapping and for the diagnostic cross-spectrum.

Use the fixed exponent `p=2` mark. For the two sign partners, define

\[
F_\pm=(1\pm t)(1\pm s)^2,\qquad
\mu_\pm=\langle F_\pm\rangle,\qquad
\delta_{m,\pm}=F_\pm/\mu_\pm-1.
\]

Each sign partner is centered by its **own** mean. Average their auto-spectra. The angular marked correlation is `M(theta) = 1 + xi_delta_m,delta_m(theta)`. Its excess over one is the Legendre transform of the marked auto-spectrum, using the output multipole `L`. Here `ell` labels input temperature modes, while `L` labels modes of the product fields and the potential.

## Measured spectra and residuals

The main harmonic descriptor is the sign-averaged marked auto-spectrum,

\[
C_L^{\rm mark}={1\over2}\left(C_L^{\delta_{m,+}\delta_{m,+}}+
C_L^{\delta_{m,-}\delta_{m,-}}\right).
\]

For the explicit connected four-point projection, form `q = t s - <t s>` and subtract its disconnected Wick contribution,

\[
T_{q,L}=C_L^{qq}-W_L^{qq},\qquad
W^{qq}(\theta)=\xi_{tt}(\theta)\xi_{ss}(\theta)+\xi_{ts}(\theta)^2.
\]

The Wick term is projected to `L` by a full-angle Legendre integral. In the sign-averaged marked auto-spectrum, the coefficient of `C_L^qq` is the realization-dependent `K_q=2(mu_+^{-2}+mu_-^{-2})`. In the present `ell >= 100` realization, `K_q` is approximately 0.271. Plot `K_q T_q` as the explicitly identified marked component, alongside the full marked result and their difference. The full mark also contains other products and normalization terms.

For each map, one phase-randomized control retains its measured temperature power spectrum. Subtract that control's **own marked and connected spectra** from the map's spectra. Then subtract the unlensed residual from the lensed residual. This yields the paired marked and connected responses shown in the pilot figures. The current control count is one per map, so these signed one-realization responses have no estimated significance.

The diagnostic `C_L^{q phi}` crosses the quadratic temperature field with the known input potential. Its paired lensed-minus-unlensed value checks where `q` tracks the injected lensing pattern. It probes a three-field response. The connected `q` auto-spectrum above is the four-field projection needed for the trispectrum study.

## Output scales and what the pilot showed

The saved spectra cover `40 <= L < 3000`. The current 88 plotting bins have width 10 over `40 <= L < 200`, width 25 over `200 <= L < 1000`, and width 50 over `1000 <= L < 3000`. Bandpowers use weights `2L+1`. Unbinned spectra are also saved. The bins of width five used in the diagnostic are a rebinning of those same unbinned spectra.

At `theta=2 arcmin`, the paired connected-`q` correlation is `-2.61e-3` in this realization. The exact multipole sum over `100 <= L < 600` contributes `-2.57e-3`, so this output range is a useful *pilot inspection range*. At `theta=30 arcmin`, different output ranges cancel and the paired connected correlation is `-3.14e-5`. Using the current band means in place of every unbinned `C_L` changes the 2-arcmin correlation by about 0.06 percent. The angle-space Wick calculation and harmonic calculation recover the same residual spectrum to `2.7e-14` in absolute spectrum units.

The present `ell >= 100` filter is the **working baseline chosen for this exercise**. The pilot does not identify a statistically best filter, smoothing width, angular range, or output `L` range. Adding `100 <= ell < 600` changes about 81 percent of the retained input temperature variance in the CAMB model and changes the mark normalization strongly. The `ell >= 600` comparison gave a more coherent `q`–`phi` cross-spectrum over a wider output range, while the one-map marked phase residual changed sign between the two filters. That contrast cannot rank their trispectrum sensitivity with one map and one phase control.

## Next test

Keep the `ell >= 100`, `sigma_s=2 arcmin`, `p=2` configuration fixed for a small paired multi-seed pilot. Use several phase controls per map, then estimate the null distribution and covariance of the signed marked residual and `K_q T_q`. Inspect the unbinned or narrow-bin output around `100 <= L < 600`, as well as the exact angular sums at a small predefined set of separations such as 2, 5, and 10 arcmin. Test resolution dependence of those small-angle results before interpreting their amplitude. Compare the `ell >= 600` temperature cut on the same seeds as a separate filter check. Choose a primary range only from ensemble consistency and covariance, not the amplitude of this single realization.
