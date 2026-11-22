# Nonlinear marked correlations as joint compressed probes of the bispectrum and trispectrum

## Aim

Marked statistics turn a nonlinear function of the density field into a two-point measurement.  This makes them a compact way to retain selected information from higher-order correlation functions without measuring the full bispectrum or trispectrum.

The central question is whether a nonlinear marked *auto*-correlation captures information that is complementary to the linear-mark cross-spectrum studied by Ebina, White & Chaussidon (2026), after the ordinary power spectrum and the full covariance are accounted for.

This note uses a deliberately controlled comparison: a non-Gaussian (NG) field and a Gaussianized field with the same two-point function.  Any difference in the marked statistic then comes from higher-order information.

## Fields, smoothing, and normalization

Let \(\delta(\mathbf{x})\) be the density contrast and define an environmental field smoothed on scale \(R\),

\[
s(\mathbf{x}) \equiv \delta_R(\mathbf{x})
= \int d^3y\,W_R(\mathbf{x}-\mathbf{y})\,\delta(\mathbf{y}),
\qquad
\delta_R(\mathbf{k})=W_R(k)\delta(\mathbf{k}).
\]

The smoothing scale is part of the estimator: it defines the environment entering the mark, suppresses sensitivity to poorly modelled UV modes, and selects a particular windowed projection of the bispectrum.

For a density-weighted marked field, define

\[
F(\mathbf{x}) \equiv [1+\delta(\mathbf{x})]m[s(\mathbf{x})],
\qquad
M(r)\equiv \frac{\langle F_1F_2\rangle}{\langle F\rangle^2}.
\]

Here \(M(r)\) is the mean-normalized weighted-pair statistic.  If one instead adopts the common marked-correlation ratio convention, divide it by the unweighted pair factor:

\[
\mathcal{M}(r)=\frac{M(r)}{1+\xi_{\delta\delta}(r)}.
\]

The distinction does not alter the higher-order interpretation, but it should be kept fixed when comparing estimators.

## Ebina's linear-mark cross-spectrum

Ebina, White & Chaussidon use a density-weighted marked field with a linear environmental mark,

\[
\frac{m(\mathbf{x})}{\bar m}=C_0+C_1s(\mathbf{x}),
\qquad
\delta_M=\frac{m}{\bar m}(1+\delta)-1.
\]

Although the mark is linear, the marked density is not:

\[
\delta_M
\supset C_1\,\delta s.
\]

Consequently, the unmarked--marked cross-correlation contains, after its two-point contribution is removed,

\[
\langle \delta_1\delta_{M,2}\rangle_{\rm beyond\ 2pt}
\propto C_1\langle\delta_1\delta_2s_2\rangle.
\]

This is a collapsed, smoothed three-point function,

\[
\zeta_{\delta,\delta s}(r)
\equiv
\langle\delta(\mathbf{x})\,\delta(\mathbf{x+r})\,\delta_R(\mathbf{x+r})\rangle.
\]

In Fourier space it is an integrated bispectrum.  Schematically,

\[
\zeta_{\delta,\delta s}(K)
\propto
\int \frac{d^3q}{(2\pi)^3}
W_R(|\mathbf{K}+\mathbf q|)
B(\mathbf K,\mathbf q,-\mathbf K-\mathbf q).
\]

At \(K\ll q\), this is a squeezed-bispectrum projection.  It is a particularly clean beyond-power-spectrum baseline because it isolates one well-defined cubic operator.

## Nonlinear environmental mark

Consider instead

\[
m(s)=(1+s)^p
=1+p_1s+p_2s^2+p_3s^3+\mathcal O(\delta^4),
\]

with

\[
p_1=p,
\qquad
p_2=\frac{p(p-1)}{2},
\qquad
p_3=\frac{p(p-1)(p-2)}{6}.
\]

Expanding the weighted density by total field order gives

\[
F=1+L+Q+C+\mathcal O(\delta^4),
\]

\[
L=\delta+ps,
\qquad
Q=p\,\delta s+p_2s^2,
\qquad
C=p_2\,\delta s^2+p_3s^3.
\]

The normalization by \(\langle F\rangle^2\) removes the purely local variance and skewness pieces.  Through third order,

\[
M(r)=1+\langle L_1L_2\rangle
+\langle L_1Q_2\rangle+\langle Q_1L_2\rangle
+\mathcal O(\delta^4).
\]

The two-point contribution is

\[
\langle L_1L_2\rangle
=\xi_{\delta\delta}(r)
+2p\,\xi_{\delta s}(r)
+p^2\xi_{ss}(r).
\]

Define collapsed three-point functions by

\[
\zeta_{a,bc}(r)
\equiv
\langle a(\mathbf{x})b(\mathbf{x+r})c(\mathbf{x+r})\rangle.
\]

Then the full cubic contribution of this marked auto-correlation is

\[
\begin{aligned}
M(r)={}&1+\xi_{\delta\delta}
+2p\,\xi_{\delta s}
+p^2\xi_{ss}\\
&+2\left[
p\,\zeta_{\delta,\delta s}
+p_2\,\zeta_{\delta,s^2}
+p^2\,\zeta_{s,\delta s}
+pp_2\,\zeta_{s,s^2}
\right]
+\mathcal O(\delta^4).
\end{aligned}
\]

The Ebina projection \(\zeta_{\delta,\delta s}\) is therefore contained in the nonlinear marked auto-correlation, but it is accompanied by three additional, differently windowed bispectrum projections.

As a consistency check, if \(R\to0\), then \(s=\delta\) and

\[
M(r)=1+(p+1)^2\xi(r)+p(p+1)^2\zeta_{21}(r)+\mathcal O(\delta^4),
\]

where \(\zeta_{21}(r)=\langle\delta^2(\mathbf{x})\delta(\mathbf{x+r})\rangle\).

## Matched-power-spectrum difference

Let the NG and Gaussianized fields have the same \(P(k)\), and apply the same smoothing filter \(W_R\).  Then \(\xi_{\delta\delta}\), \(\xi_{\delta s}\), and \(\xi_{ss}\) are identical.  Since all cubic correlators vanish for the Gaussian field,

\[
\begin{aligned}
\Delta M(r)
&\equiv M_{\rm NG}(r)-M_{\rm G}(r)\\
&=2\left[
p\,\zeta_{\delta,\delta s}^{\rm NG}
+p_2\,\zeta_{\delta,s^2}^{\rm NG}
+p^2\,\zeta_{s,\delta s}^{\rm NG}
+pp_2\,\zeta_{s,s^2}^{\rm NG}
\right]
+\mathcal O(\delta^4).
\end{aligned}
\]

Thus \(\Delta M\) is a set of compressed, windowed bispectrum projections at leading NG order.  It is not the full squeezed bispectrum: each term integrates over the short-scale magnitude and angular dependence of \(B(K,q,\mu)\).

## What is gained relative to the Ebina cross-statistic?

At fixed \(p\) and \(R\), the nonlinear marked auto-correlation is still one summary statistic.  It mixes several projections; this does **not** automatically imply more independent information than the Ebina cross-spectrum.

Nevertheless, it has two clear routes to complementarity.

1. **Different bispectrum windows.**  The Ebina statistic measures one inner product of \(B\) with a kernel \(\mathcal W_{\rm E}\).  The nonlinear auto-statistic uses other kernels involving one, two, or three smoothing windows.  A scale-dependent bispectrum response can be orthogonal to \(\mathcal W_{\rm E}\), giving a vanishing Ebina signal, while remaining non-zero under one of the nonlinear weights.

2. **Trispectrum sensitivity.**  At fourth order, the marked auto-correlation contains terms such as

   \[
   \langle(\delta s)_1(\delta s)_2\rangle_c,
   \qquad
   \langle s_1^2s_2^2\rangle_c,
   \]

   which probe connected trispectrum projections.  The linear unmarked--marked cross-spectrum is a cubic statistic and has no corresponding leading trispectrum term.  A parity-symmetric non-Gaussian field with \(B=0\) but \(T\neq0\) provides a decisive toy example: after matching \(P(k)\), the Ebina cross residual vanishes while the marked auto-correlation can remain non-zero.

The second point is principally an advantage of including a marked *auto*-correlation; nonlinear marks supply additional composite operators and therefore additional trispectrum windows.

## Reproducible trispectrum-only tests

The repository contains two controlled tests of the statement above.  They
are not intended to model a late-time density field or a calibrated CMB
cosmology; their purpose is to isolate the estimator response when the
bispectrum is absent and a connected four-point function is present.

### Periodic-grid cubic field

`lss/trispectrum_marked_comparison.py` first draws a Gaussian periodic field
\(g\), then applies the odd local transform

\[
\delta=g+A\left(g^3-3\langle g^2\rangle g\right).
\]

The probability law is invariant under \(\delta\to-\delta\), hence its
ensemble odd moments, including the bispectrum, vanish.  The cubic transform
nevertheless generates a connected trispectrum.  To eliminate finite-volume
odd-order leakage exactly, every statistic is averaged over \(\delta\) and
\(-\delta\).  A phase-randomized copy of \(\delta\), paired in the same way,
is the matched-power control.

The nonlinear statistic uses the \((1+\delta)(1+\delta_R)^p\) mark already
defined above.  The comparison statistic uses Ebina's linear mark,

\[
\delta_{\rm E}=\frac{(1+\delta)(1+\delta_R)}
{\langle(1+\delta)(1+\delta_R)\rangle}-1,
\qquad
X_{\rm E}(r)=\langle\delta_1\delta_{{\rm E},2}\rangle.
\]

The raw \(X_{\rm E}\) is non-zero because it contains ordinary two-point
terms.  Its difference from the matched-power control is zero in the paired
test, whereas the nonlinear marked-auto residual is non-zero.  For the
default \(64^3\) realization, the maximum nonlinear residual is
\(6.28\times10^{-3}\), while the Ebina residual is at floating-point zero
(\(\lesssim5\times10^{-17}\)).  This is the direct numerical realization of
\(B=0,\ T\neq0\) sensitivity.

Both the periodic-grid and CMB pipelines also report an explicit connected
four-point estimator.  Define the centered quadratic operator

\[
Q=\delta s-\langle\delta s\rangle,
\qquad C_Q(r)=\langle Q_1Q_2\rangle.
\]

Its Wick decomposition gives an exact collapsed trispectrum projection,

\[
T_Q(r)=C_Q(r)-\xi_{\delta\delta}(r)\xi_{ss}(r)
-\xi_{\delta s}(r)\xi_{s\delta}(r).
\]

Since the phase-randomized control has the same two-point functions, the
implementation checks \(\Delta C_Q=\Delta T_Q\).  The default grid run gives
\(\max_r|\Delta C_Q-\Delta T_Q|=5.4\times10^{-20}\), directly verifying that
this residual is the connected trispectrum projection rather than a
two-point or Wick-disconnected contribution.

The accompanying `trispectrum_decomposition.png` compares the full marked
residual with the contribution \(p^2\Delta T_Q\) implied by the
\(p\,\delta s\) term in the quadratic marked operator.  Their difference is
not an error term: it contains the other trispectrum projections in
\(Q_{\rm mark}Q_{\rm mark}\) and \(LC+CL\), as well as higher even-order
terms beyond the weak-non-Gaussian expansion.

This gives a stronger validation than observing two non-zero statistics.  The
direct estimator proves algebraically and numerically that one known component
of the marked residual is a connected trispectrum contraction: in the default
run, \(\Delta C_Q=\Delta T_Q\) to \(5.4\times10^{-20}\).  The decomposition
then predicts that this component enters the full marked statistic with the
fixed coefficient \(p^2\), without fitting an amplitude.  The observed
\(\Delta M\) follows this contribution on large scales and has a remaining
response of at most \(2.80\times10^{-3}\), consistent with the other
explicitly identified composite-operator contractions and, at finite
amplitude, higher even-order terms.  Thus the test
connects the nonlinear marked response to the trispectrum structurally, not
merely by a shared non-zero signal.

### Full-sky CMB lensing field

`archive/cmb_pilots/scripts/cmb_lensing_trispectrum_demo.py` generates a HEALPix temperature map using
an independent Gaussian primary temperature \(T\) and lensing potential
\(\phi\), with \(C_\ell^{T\phi}=0\).  Lensing remaps the primary field and
therefore generates a connected CMB trispectrum.  The independence and
\(T\to-T\) symmetry ensure a vanishing ensemble temperature bispectrum.  The
default spectra are deliberately simple toy spectra; replacing them with CAMB
\(C_\ell^{TT}\) and \(C_\ell^{\phi\phi}\), while retaining
\(C_\ell^{T\phi}=0\), gives the corresponding clean cosmological test.

`archive/cmb_pilots/scripts/cmb_marked_statistics.py` implements the field-level angular analogue.  Let
\(t=T/\sigma_T\), and let \(s\) be its Gaussian-harmonically smoothed map.  It
uses

\[
m_{\rm NL}(s)=1+a s+b\left(s^2-\langle s^2\rangle\right),
\qquad \widetilde F_{\rm NL}=t\,m_{\rm NL}(s)-\langle t\,m_{\rm NL}(s)\rangle,
\]

and measures \(C_{\rm NL}(\theta)=\langle\widetilde F_{{\rm NL},1}\widetilde F_{{\rm NL},2}\rangle\).
The CMB linear-mark analogue is

\[
X_{\rm E}(\theta)=
\left\langle t_1\,[t(1+a s)-\langle t(1+a s)\rangle]_2\right\rangle.
\]

The same direct diagnostic is evaluated on the sphere with
\(Q=t s-\langle t s\rangle\).  Replacing every \(\xi\) in the preceding
Wick subtraction by its angular correlation gives \(T_Q(\theta)\).  In the
supplied test, \(\Delta C_Q(\theta)\) and \(\Delta T_Q(\theta)\) agree to
better than \(2\times10^{-8}\); the small mismatch is spherical-transform and
pixelization accuracy, rather than a physical contribution.

The control is made by randomizing the phases of the spherical-harmonic
coefficients \(a_{\ell m}\), retaining every \(|a_{\ell m}|\) and thus the
angular power spectrum \(C_\ell\) exactly.  Sign-pairing \(t\) and \(-t\)
again cancels odd-order contributions in the reported estimators.  The code
therefore compares the nonlinear and linear-mark angular correlations on a
common \(\theta\) axis, while the phase-randomized difference of
\(C_{\rm NL}\) is the lensing-trispectrum response.  The supplied full-sky
test has matched \(C_\ell\) values to \(4\times10^{-16}\).

For harmonic-space interpretation, the same code additionally writes
`cmb_trispectrum_decomposition_ell.png`.  It Legendre-transforms the
configuration-space Wick-subtracted \(T_Q(\theta)\) and compares
\(\Delta C_\ell^{\rm NL}\) with the predicted \(a^2\Delta C_{\ell,T_Q}\)
component on an \(\ell\) axis.  This is the appropriate decomposition plot
for the CMB map; the separate periodic-grid decomposition is correctly a
function of \(r\).

For a masked sky, anisotropic noise, or beam and transfer-function treatment,
this full-sky estimator should be replaced by a pseudo-\(C_\ell\) pipeline
such as NaMaster, with identically processed simulations providing the
mean-field and disconnected-term corrections.

## A falsifiable study programme

The study should test whether the nonlinear estimator adds independent constraining power, rather than merely displaying a different non-zero signal.

### 1. Analytic operator basis

Derive the marked cross- and auto-correlations through fourth order, including the connected and disconnected four-point contributions.  Identify which terms are removed when comparing matched-\(P(k)\) fields.

### 2. Controlled field tests

Use three field classes:

- original late-time density fields and phase-randomized fields with matched \(P(k)\), to isolate ordinary gravitational NG;
- synthetic fields with a controlled scale-dependent bispectrum, to construct cancellations in the Ebina window but not in nonlinear windows;
- parity-symmetric non-Gaussian fields with \(B=0\) and \(T\neq0\), to isolate trispectrum sensitivity.

### 3. Cosmological information comparison

Measure the joint data vector

\[
\mathbf d=
\left\{
P_{\delta\delta},
\ P_{\delta\,\delta_M}^{\rm linear}(R),
\ P_{\delta_M\delta_M}^{\rm linear}(R),
\ P_{\delta_M\delta_M}^{\rm nonlinear}(p,R)
\right\}.
\]

Vary \(R\) and \(p\), and use the full joint covariance.  Fisher forecasts or simulation-based inference can then determine whether the nonlinear auto-statistic improves constraints beyond the ordinary power spectrum plus the Ebina-style cross-spectrum.

## Falsifiable outcomes and interpretation

- If the nonlinear auto-mark is fully predicted by \(P\) plus the linear cross-spectrum, it is a redundant compression in the tested regime.
- If it improves parameter constraints only for particular \(p\) or \(R\), the associated kernels identify which environmental scales carry the additional information.
- If it responds in \(B=0,\ T\neq0\) tests where the linear cross-spectrum does not, it demonstrates direct sensitivity to higher-than-bispectrum information.

Either outcome is informative.  The interesting result is not simply that nonlinear transformations are non-Gaussian, but whether a controlled family of marked auto-correlations gives a compact, robust, and interpretable extension of existing bispectrum-compressed marked statistics.

## Relevant literature

- M. White, *A marked correlation function for constraining modified gravity models*, JCAP **11** (2016) 057, [arXiv:1609.08632](https://arxiv.org/abs/1609.08632).
- O. H. E. Philcox, E. Massara & D. N. Spergel, *What does the marked power spectrum measure? Insights from perturbation theory*, Phys. Rev. D **102** (2020) 043516, [arXiv:2006.10055](https://arxiv.org/abs/2006.10055).
- H. Ebina, M. White & E. Chaussidon, *The marked power spectrum as a practical bispectrum measure for galaxy redshift surveys*, JCAP **08** (2026) 003, [arXiv:2603.03581](https://arxiv.org/abs/2603.03581).
- J. A. Cowell, D. Alonso & J. Liu, *Optimizing marked power spectra for cosmology*, MNRAS **535** (2024) 3129, [DOI: 10.1093/mnras/stae2492](https://doi.org/10.1093/mnras/stae2492).
