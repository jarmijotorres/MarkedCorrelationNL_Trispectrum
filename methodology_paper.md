# Methodology: testing trispectrum sensitivity of nonlinear marked correlations

## Scope

Test whether a fixed nonlinear marked auto-correlation responds to connected four-point structure when the ensemble bispectrum vanishes, while a linear-mark cross-correlation has no corresponding higher-order response. Use two settings: a periodic 3D field and a full-sky lensed CMB temperature field. The claim is trispectrum sensitivity. Comparison with a CMB lensing quadratic estimator, optimization of marks, and incremental cosmological information are reserved for later work.

The existing plots and the ten-pair Nside=128 CMB benchmark are algorithm checks. Paper-level conclusions require an input-spectrum and scale-sensitivity check before a large ensemble, followed by independent ensembles, null distributions, uncertainty bands, and numerical convergence checks.

## Fixed estimators

### Periodic 3D field

For mean-zero \(\delta(\mathbf x)\), smooth \(s(\mathbf k)=W_R(k)\delta(\mathbf k)\), with \(W_R(k)=\exp(-k^2R^2/2)\). Define

\[
F_p=(1+\delta)(1+s)^p,\qquad f_p=F_p/\langle F_p\rangle-1,\qquad
M_p(r)=1+\langle f_{p,1}f_{p,2}\rangle .
\]

Fix \(p=2\) and \(R=10\) for the initial toy box. This is a mean-normalized weighted-pair statistic, not the ratio \(M_p/(1+\xi_{\delta\delta})\). The linear-mark comparison is \(X_E(r)=\langle\delta_1 f_{1,2}\rangle\), with \(f_1\) using the same normalization.

For the direct connected four-point projection, define \(q=\delta s-\langle\delta s\rangle\):

\[
C_q(r)=\langle q_1q_2\rangle,\qquad
T_q(r)=C_q(r)-\xi_{\delta\delta}(r)\xi_{ss}(r)
-\xi_{\delta s}(r)\xi_{s\delta}(r).
\]

Evaluate the Wick products at each grid-separation vector before radial averaging, then bin all terms identically. The current implementation multiplies separately binned correlations, an approximation to the bin-averaged Wick term.

### Full-sky CMB temperature field

Define a dimensionless, variance-normalized temperature contrast \(\delta_T=(T-\langle T\rangle)/\sigma_T\), separately for each filtered map. This \(\delta_T\) is not the physical CMB contrast \((T-\bar T)/\bar T\); the normalization is an explicit toy choice. Smooth the temperature contrast as \(s_{\ell m}=W_\ell(\delta_T)_{\ell m}\), where \(W_\ell=\exp[-\ell(\ell+1)\sigma_s^2/2]\). The historical low-resolution \(N_{\rm side}=128\) toy benchmark used \(\sigma_s=2^\circ\); the current \(N_{\rm side}=2048\) CMB preflight uses \(\sigma_s=2\) arcmin. This is the Gaussian standard deviation, corresponding to FWHM \(4.71\) arcmin. The lensing potential \(\phi\) is not smoothed in this step. Use the same \(p=2\) mark as in the 3D test:

\[
m(s)=(1+s)^2,\qquad
F=(1+\delta_T)m(s),\qquad
\delta_m^{\rm CMB}=\frac{F}{\langle F\rangle}-1.
\]

Measure \(M_{\rm CMB}(\theta)=1+\langle\delta_{m,1}^{\rm CMB}\delta_{m,2}^{\rm CMB}\rangle\) and its \(C_\ell^{\delta_m\delta_m}\). Since CMB temperature is not a density field, \(1+\delta_T\) is an algebraic analogue of the density prefactor, not a physical pair-count weight. For \(q=\delta_T s-\langle\delta_T s\rangle\), compute

\[
T_q(\theta)=C_q(\theta)-\xi_{\delta_T\delta_T}(\theta)\xi_{ss}(\theta)
-\xi_{\delta_Ts}(\theta)\xi_{s\delta_T}(\theta).
\]

For sign pairing, write \(F_\pm=(1\pm\delta_T)(1\pm s)^2=1+y\pm x\), with \(x=\delta_T+2s+\delta_Ts^2\), \(y=2\delta_Ts+s^2\), and \(\mu_\pm=1+\langle y\rangle\pm\langle x\rangle\). Normalize the two sign partners separately. In the sign-averaged marked auto-spectrum, the \(C_\ell^{qq}\) term has the exactly known, realization-dependent coefficient \(K_q=2(\mu_+^{-2}+\mu_-^{-2})\). The displayed connected piece is \(K_qT_q\) for each map and each phase control before subtraction. It is not a fitted constant and need not equal four at unit variance. Legendre-transform the Wick projection for the \(\ell\)-space decomposition. Record \(N_{\rm side}\), \(\ell_{\max}\), pixel window, smoothing, harmonic iterations, and plotted scale cuts. The baseline is full sky and noiseless.

## Desired ensemble and null test

Use independent simulation seeds as the units of replication. Pair each injection with a null generated from the same Gaussian seed. Fix the number of seeds and analysis bins before inspecting ensemble results; increase the sample only to stabilize the mean uncertainties. Archive the configuration and seeds.

1. **3D null and injection.** Draw Gaussian \(g\). The null has \(A=0\); the injection uses \(\delta_A=g+A(g^3-3\langle g^2\rangle g)\), initially with the current fixed \(A=2\). The sign-symmetric ensemble has vanishing odd moments and a connected four-point response for nonzero \(A\). Sign-pair every estimator on \(\delta\) and \(-\delta\) to cancel finite-realization odd terms. These partners constitute one measurement, not two independent volumes. Verify that the density-style mark and its normalization are well defined.
2. **CMB null and injection.** Draw independent Gaussian primary \(T\) and Gaussian lensing potential \(\phi\) with \(C_\ell^{T\phi}=0\). Compare the unlensed primary with the same primary lensed by lenspyx, and sign-pair the estimators. Since lensing changes the realized \(C_\ell\), the lensed-minus-unlensed difference is not by itself a pure trispectrum difference.
3. **Matched two-point controls.** For each injected 3D or CMB field, generate multiple independent phase-randomized surrogates preserving its exact Fourier or harmonic amplitudes and real-field conditions. Verify \(P(k)\) or \(C_\ell\) matching. Use surrogate scatter in uncertainties. A finite phase surrogate is not guaranteed to be Gaussian, and surrogates of one injection are not independent injected realizations.
4. **Independent Gaussian checks.** Process the \(A=0\) and unlensed CMB ensembles with identical smoothing, centering, sign pairing, binning, and transforms. Verify that their connected projection and matched-spectrum marked residuals agree with their empirical null distributions. Keep CMB secondary effects capable of generating a bispectrum out of this controlled test.

The two settings use distinct coordinates: \(r\) in the box, \(\theta\) for angular correlations, and \(\ell\) for CMB harmonic spectra. Do not compare their numerical amplitudes directly.

## Mechanism tests

**Bispectrum null.** Measure the collapsed cubic projection entering the Ebina-style cross-statistic on the *unpaired* ensemble; its mean should be consistent with zero. Sign pairing also cancels the odd contribution to the measured cross residual exactly. The unpaired check validates the simulator's intended symmetry.

**Connected four-point response.** Show that the injected ensemble has nonzero \(T_q\) relative to its matched-spectrum controls, while the Gaussian null is consistent with zero. Check \(\Delta C_q=\Delta T_q\) after correcting the 3D binning. The identity follows from matched disconnected Wick terms; its tiny single-realization error is a numerical check, not detection significance.

**Predicted marked component.** Derive the centered, mean-normalized mark through fourth order, including normalization terms. The \(p\,\delta s\) term in the 3D quadratic operator has nominal coefficient \(p^2=4\) before normalization. The corresponding CMB \(q=\delta_Ts-\langle\delta_Ts\rangle\) auto-term has the exact normalization coefficient \(K_q=2(\mu_+^{-2}+\mu_-^{-2})\) in the sign-paired spectrum. Plot signed ensemble means of the full marked residual, \(K_qT_q\) after matched-control subtraction, and the remainder, with seed-based uncertainty bands. Do not fit a multiplier. The remainder includes other connected four-point contractions, normalization effects, and higher even orders.

**Convergence.** Repeat the fixed-parameter test at higher mesh or \(N_{\rm side}\)/\(\ell_{\max}\) to establish stable scales. A small-injection check can verify perturbative behavior if needed; no mark-parameter scan is required.

**CMB scale-sensitivity gate.** In each case, \(s\) is the smoothed *lensed or unlensed temperature map*, respectively; it is never a smoothed lensing-potential map. Distinguish lensing-potential multipoles \(L\) from temperature multipoles \(\ell\). The Nside=128, \(\ell_{\max}=128\), \(\sigma_s=2^\circ\) benchmark has not established sensitivity to realistic lensing: it only includes lensing modes through \(L=128\), truncates the temperature field at \(\ell=128\), and uses toy power-law spectra. Its roughly 27.5-arcmin pixels are larger than the roughly 2-arcmin toy deflection, but subpixel deflection alone does not prove insensitivity. The 2-degree window keeps mainly \(\ell\lesssim40\) temperature modes in \(s\); \(q=\delta_Ts\) can still couple them to finer temperature modes in \(\delta_T\). A lensing mode around \(L=100\) can contribute, while \(L\sim200\) is absent from the present output range. Do not treat the \(\ell\sim2000\) reach of a standard quadratic reconstruction as a mandatory cutoff for this different descriptor. Before enlarging the ensemble, use physical CMB spectra and compare the lensing-induced four-point response with the unlensed null over a small resolution/\(\ell_{\max}\) ladder that includes \(L\sim100\)–200. Fix the smoothing window and analysis range only after this response check; a small sky patch can provide a cheaper preflight than a full-sky thousand-map run.

**ACT DR6 reference choices.** ACT's lensing reconstruction uses input CMB modes \(600<\ell<3000\), a baseline lensing-spectrum range \(40<L<763\), and 1-arcmin downgraded map pixels. Its first four reported lensing bins are \(L=40\)–66, 67–101, 102–145, and 146–199. The collaboration tested several alternative input-\(\ell\) cuts, including 300–2000, and generated lensed Gaussian-CMB simulations from fiducial theory spectra. These choices motivate checking our response around \(L\sim100\)–200, using physical spectra, and varying the temperature-mode cut for convergence. They do not prescribe \(\sigma_s=2^\circ\) for our mark: applying that Gaussian to an ACT-filtered map would suppress all retained \(\ell\ge600\) modes (\(W_{600}\approx4\times10^{-96}\)). An ACT-like application must define a new temperature filter for \(s\), or supply independent large-scale temperature information. Also, the multipole of our \(C_\ell^{qq}\) projection is not automatically an unbiased measurement of lensing \(L\); calibrate its response with input-\(\phi\) simulations before equating the axes. Sources: [ACT DR6 map and methods](https://ntrs.nasa.gov/api/citations/20240002550/downloads/Madhavacheril_ACT%20DR6%20Lensing%20Cos%20Parameters%20final.pdf), [ACT DR6 lensing spectrum](https://lss.fnal.gov/archive/2023/pub/fermilab-pub-23-237-ppd.pdf), [ACT simulation release](https://lambda.gsfc.nasa.gov/product/act/actadv_dr6_lensing_maps_info.html).

**Current CMB preflight choice.** Use an input-temperature passband \(100\le\ell\le3000\) and a Gaussian temperature mark window with \(\sigma_s=2\) arcmin (FWHM about 4.7 arcmin), not 2 degrees. This window has half-amplitude at \(\ell\approx2023\): \(W_{600}=0.94\), \(W_{1500}=0.68\), \(W_{2000}=0.51\), and \(W_{3000}=0.22\). Retaining \(100\le\ell<600\) supplies large-scale temperature information; the earlier ACT-like \(600\le\ell\le3000\) run remains a controlled filter comparison. Measure the marked and \(q\) output spectra beyond \(L=200\), and compare their angular correlations at arcminute separations. Use physical fiducial CMB and lensing spectra. For a full-sky pilot at \(\ell_{\max}=3000\), choose HEALPix \(N_{\rm side}=2048\) (about 1.72-arcmin pixels) as the starting resolution. A real scalar-map product can contain multipoles up to roughly 6000, close to the nominal HEALPix transform limit \(3N_{\rm side}-1=6143\), so validate the marked and \(q\) results against a finer resolution, such as \(N_{\rm side}=4096\), on a few matched realizations before running an ensemble. Because full-sky maps at these resolutions are expensive, a 1-arcmin flat-sky patch is a practical first preflight. This is a declared filter for the response test, not an optimized mark or a claim that our output multipole equals lensing \(L\).

The first full-sky, noiseless `Nside=2048` realization using Planck-like CAMB spectra is saved in `archive/cmb_pilots/run_cmb_actlike_single_2048/`. Its `q` field correlates with the known input lensing potential at about -0.06 to -0.08 across the four low output bins, while its paired unlensed field is nearer zero. This is a single-map response diagnostic; the marked auto-spectrum response and uncertainty remain unresolved. The full paired run with one phase control per map took about two minutes, making a large full-sky ensemble costly without further optimization or a patch-based approach.

An extended-output preflight using the same paired seed is saved in `archive/cmb_pilots/run_cmb_actlike_high_ell_2048/`. The input temperature range remains \(600\le\ell\le3000\); the measured marked and quadratic outputs now cover \(40\le L<3000\). The angular correlation from only \(L<200\) omitted most of the one-map paired marked response below about 15 arcmin: at 2 arcmin the band-limited values are \(4.19\times10^{-4}\) for \(L<200\) and \(62.8\times10^{-4}\) for \(L<3000\). The normalized connected-\(q\) component is \(1.79\times10^{-4}\) versus \(10.3\times10^{-4}\). These are single-realization, phase-control-subtracted values, not detection significances or an isolated trispectrum measurement. The cross-spectrum \(C_L^{q\phi}\) confirms a lensing response over a wider \(L\) range but is a three-field diagnostic. For the paper, retain output scales beyond \(L=200\) and small angular separations, and determine their covariance and resolution dependence with an ensemble. The product fields can contain modes above the current output cap \(L=3000\), so test convergence before calling the angular correlation complete.

The matching \(100\le\ell\le3000\) preflight is saved in `archive/cmb_pilots/run_cmb_actlike_lmin100_2048/`. Modes \(100\le\ell<600\) contribute about 81% of the input temperature variance in the CAMB spectrum, so adding them changes the unit-variance temperature contrast and the mark weighting substantially. For this seed, the paired \(q\)–\(\phi\) difference is coherently negative through roughly \(L=200\) and loses that sign coherence above a few hundred, while the marked and connected-\(q\) phase-control residuals fluctuate strongly. The cross-spectrum is a useful response locator in output \(L\), not an estimator of the marked trispectrum signal-to-noise. An ensemble must establish the useful \(L\) range for the marked statistic separately.

The first marked-spectrum filter screen is in run_cmb_marked_pilot_2048/ and run_cmb_marked_pilot_lmin600_2048/. It compares the same three paired skies with three phase controls per map for input cuts \(\ell\ge100\) and \(\ell\ge600\), keeping \(N_{\rm side}=2048\), \(\sigma_s=2'\), and \(p=2\) fixed. The signed paired marked bandpower over output \(100\le L<600\) changes sign across the \(\ell\ge100\) seeds, while all three \(\ell\ge600\) values are positive at \(6.19\), \(7.09\), and \(6.97\times10^{-8}\). The corresponding marked angular curves are also positive at \(2'\). The \(\ell\ge600\) input is the leading candidate for the marked observable. This three-seed screen has no established control covariance or resolution convergence, so it does not yet support a detection claim. The direct connected-\(q\) term remains a decomposition check for the marked result.

**Angular and harmonic consistency.** The repeat in `run_cmb_angular_consistency_2048/` saved the raw temperature and \(q\) spectra for each map and phase control. We rebuilt \(T_q(\theta)=\xi_{qq}(\theta)-\xi_{TT}(\theta)\xi_{ss}(\theta)-\xi_{Ts}(\theta)^2\), with \(T\) denoting the filtered temperature contrast in these correlations. We applied the same phase and lensed-minus-unlensed subtractions, then integrated against Legendre polynomials over the full angular range. The reconstructed residual spectra agree with the harmonic calculation to \(2.7\times10^{-14}\) in absolute spectrum units over \(40\le L<3000\). At \(30'\), the contributions from \(40\le L<200\) and \(200\le L<3000\) nearly cancel in the paired correlation, giving \(-3.14\times10^{-5}\). At \(2'\), both bands contribute with the same sign and give \(-2.61\times10^{-3}\). This check establishes consistency for one realization. Use the same signed output bands and controls for comparisons, carry their covariance into any angle-space summary, and check resolution at the smallest separations.

Before a 100–1000-pair CMB benchmark, stabilize the matched-power control and make the runner checkpoint each seed. In one fixed lensed/unlensed pair, switching from independent to shared phase-control random seeds changed the 67–101 marked contrast from negative to positive; one surrogate per map is not a stable basis for the primary result. Plot raw marked and connected bandpowers rather than dividing them by \(C_\ell^{TT}\); that ratio is especially unsuitable when the input temperature filter suppresses the output-bin power. Integrate the paired known-input-\(\phi\) cross-check into the per-seed outputs and check the relevant output bins for pixelization convergence. Then run a small multi-seed variance pilot to choose the main ensemble size. Details are in `archive/cmb_pilots/run_cmb_actlike_single_2048/benchmark_readiness.md`.

## Decision criterion

Predefine displayed ranges in \(r\), \(\theta\), and \(\ell\). Use independent seeds for confidence intervals and the covariance of the measured descriptors. The primary result is supported when the injected ensemble has a resolved connected four-point response and a marked response with the predicted fixed-coefficient component, while the Gaussian null and Ebina-style higher-order cross response are consistent with zero. Report uncertainty and effect size by scale, not only the maximum bin value.

This study establishes trispectrum sensitivity of the specified marks. It does not claim optimality, a full trispectrum measurement, or independent information beyond a CMB lensing quadratic estimator.

## Deliverables and work order

Archive code revision, fixed inputs, spectra, seeds, per-seed measurements, matched-power checks, null distributions, uncertainty bands, and convergence tests. Label current figures in 'run_trispectrum_64/' as single-realization demonstrations and the earlier 'archive/cmb_pilots/run_cmb_marked_64/' polynomial-mark figures as a distinct exploratory descriptor. The corrected p=2 pilot is in 'archive/cmb_pilots/run_cmb_benchmark_128_p2/'.

1. Correct the 3D Wick subtraction at grid separation and verify it against Gaussian nulls.
2. Validate CMB lensing sensitivity at realistic temperature multipoles and resolution, then freeze the smoothing and analysis range.
3. Build paired multi-seed 3D and CMB null/injection runners with per-seed matched-spectrum surrogates.
4. Validate bispectrum nulls, power matching, connected four-point response, and numerical convergence.
5. Complete the fourth-order operator derivation and plot signed ensemble decompositions with uncertainty bands.
6. Write the result around the fixed-mark, two-ensemble trispectrum sensitivity test.
