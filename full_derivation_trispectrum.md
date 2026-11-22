# Full derivation: marked correlation response to a trispectrum

## Question and conventions

We ask whether a nonlinear marked correlation responds to connected four-point information when the bispectrum vanishes. A periodic 3D field supplies a controlled algebraic test; idealized CMB lensing supplies a physical example. This derivation identifies a contribution to the response. It does not establish optimality or independent information beyond a standard CMB lensing quadratic estimator.

The density contrast \(\delta(\mathbf x)\) has zero mean. Its smoothed version is \(s(\mathbf x)=\delta_R(\mathbf x)\), with \(s(\mathbf k)=W_R(k)\delta(\mathbf k)\). We use the mark

\[
m(s)=(1+s)^p
\]

and define the **marked density contrast** by

\[
\boxed{\delta_m(\mathbf x)=
\frac{[1+\delta(\mathbf x)]m[s(\mathbf x)]}
{\langle(1+\delta)m(s)\rangle}-1.}
\]

Thus \(\langle\delta_m\rangle=0\). Our mean-normalized marked correlation function is

\[
\boxed{M(r)=1+\langle\delta_m(\mathbf x)\delta_m(\mathbf x+\mathbf r)\rangle.}
\]

This is not the alternative marked-correlation *ratio* \(M(r)/[1+\xi_{\delta\delta}(r)]\). Subscripts 1 and 2 below denote points separated by \(r\).

## Perturbative counting and the marked contrast

Both \(\delta\) and its *linearly smoothed* version \(s\) count as first order. For example, \(\delta s\) is second order and \(\delta s^2\) is third order. Define

\[
p_2=\frac{p(p-1)}{2},\qquad
p_3=\frac{p(p-1)(p-2)}{6}.
\]

Expand the mark and then its density weighting:

\[
\begin{aligned}
m(s)&=1+ps+p_2s^2+p_3s^3+O(\delta^4),\\
(1+\delta)m(s)&=1+L+Q+C+O(\delta^4),\\
L&=\delta+ps,\\
Q&=p\,\delta s+p_2s^2,\\
C&=p_2\delta s^2+p_3s^3.
\end{aligned}
\]

Here \(L,Q,C\) have total field orders one, two and three. Because \(\langle\delta\rangle=\langle s\rangle=0\), we have \(\langle L\rangle=0\). Put \(\mu_2=\langle Q\rangle\) and \(\mu_3=\langle C\rangle\). The mean in the denominator is \(1+\mu_2+\mu_3+O(\delta^4)\). Subtract this denominator from the numerator *before* expanding the quotient:

\[
\delta_m
=\frac{1+L+Q+C+\cdots}{1+\mu_2+\mu_3+\cdots}-1
=\frac{L+(Q-\mu_2)+(C-\mu_3)+\cdots}
{1+\mu_2+\mu_3+\cdots}.
\]

The inverse denominator begins \(1-\mu_2-\mu_3+\cdots\). Since \(\mu_2\) is second order, \(-\mu_2L\) is third order; \(-\mu_3L\) and \(-\mu_2Q\) first enter at fourth order. Therefore

\[
\boxed{\delta_m
=L+(Q-\mu_2)+(C-\mu_3-\mu_2L)+O(\delta^4).}
\]

This makes the role of normalization explicit: subtracting a mean would not produce the \(-\mu_2L\) term; dividing by the mean does.

For the fixed grid-test mark \(p=2\), \(m(s)=1+2s+s^2\) exactly. Then

\[
\begin{aligned}
L&=\delta+2s,&
Q&=2\delta s+s^2,&
C&=\delta s^2,\\
\mu_2&=2\langle\delta s\rangle+\langle s^2\rangle,&
\mu_3&=\langle\delta s^2\rangle.
\end{aligned}
\]

The density-weighted numerator is exactly cubic in \(\delta,s\), but the normalized \(\delta_m\) has further perturbative terms because its denominator is inverted. To compute \(M\) through fourth order, however, the displayed third-order \(\delta_m\) is sufficient: a fourth-order term in one copy of \(\delta_m\) first pairs with a first-order term in the other copy at fifth order.

## Marked correlation through fourth order, before subtraction

Let \(q=Q-\mu_2\) and \(c=C-\mu_3\). Multiplying \(\delta_{m,1}=L_1+q_1+c_1-\mu_2L_1+\cdots\) by its copy at point 2 yields

\[
\boxed{\begin{aligned}
M(r)-1={}&
\underbrace{\langle L_1L_2\rangle}_{\text{second order}}\\
&+\underbrace{\langle L_1q_2+q_1L_2\rangle}_{\text{third order}}\\
&+\underbrace{\langle q_1q_2\rangle
+\langle L_1c_2+c_1L_2\rangle
-2\mu_2\langle L_1L_2\rangle}_{\text{fourth order}}
+O(\delta^5).
\end{aligned}}
\]

The final fourth-order term is the normalization contribution. All expectations in this section refer to **one field**, before defining a difference \(\Delta M\).

For \(p=2\), write

\[
u=\langle\delta(\mathbf x)s(\mathbf x)\rangle=\xi_{\delta s}(0),\qquad
v=\langle s(\mathbf x)^2\rangle=\xi_{ss}(0).
\]

These are same-position averages: \(u=\langle\delta_1s_1\rangle=\langle\delta_2s_2\rangle\), and similarly for \(v\), by homogeneity. In particular, \(u\) is the average of the pointwise product \(\delta s\), and \(v\) is the variance of the mean-zero field \(s\). The variance of \(\delta\) is implicit in \(\xi_{\delta\delta}(0)=\langle\delta^2\rangle\); it does not appear as a separate local normalization constant.

Assuming symmetry under exchanging the two points, the second-order term is

\[
\langle L_1L_2\rangle
=\xi_{\delta\delta}(r)+4\xi_{\delta s}(r)+4\xi_{ss}(r).
\]

The third-order term is

\[
\begin{aligned}
\langle L_1q_2+q_1L_2\rangle={}&
4\langle\delta_1\delta_2s_2\rangle
+2\langle\delta_1s_2^2\rangle\\
&+8\langle s_1\delta_2s_2\rangle
+4\langle s_1s_2^2\rangle.
\end{aligned}
\]

The three pieces at fourth order expand to

\[
\begin{aligned}
\langle q_1q_2\rangle={}&
4[\langle(\delta s)_1(\delta s)_2\rangle-u^2]\\
&+4[\langle(\delta s)_1s_2^2\rangle-uv]\\
&+[\langle s_1^2s_2^2\rangle-v^2],\\[2pt]
\langle L_1c_2+c_1L_2\rangle={}&
2\langle\delta_1\delta_2s_2^2\rangle
+4\langle s_1\delta_2s_2^2\rangle,\\
-2\mu_2\langle L_1L_2\rangle={}&
-2(2u+v)[\xi_{\delta\delta}+4\xi_{\delta s}+4\xi_{ss}].
\end{aligned}
\]

The first line contains centered four-field expectations, which still include disconnected products of two-point functions. The \(\mu_3\) subtraction in \(c\) drops out of \(\langle Lc\rangle\) because \(\langle L\rangle=0\). Terms of order five and above have been omitted.

## Defining the matched-power difference

Let an injected field and a phase-randomized control have identical Fourier-mode amplitudes. With the same linear smoothing, they have identical \(\xi_{\delta\delta}(r)\), \(\xi_{\delta s}(r)\), \(\xi_{ss}(r)\), and hence identical \(u,v,\mu_2\). Define

\[
\Delta M(r)=M_{\rm injected}(r)-M_{\rm control}(r).
\]

This is an engineered counterfactual, not the difference between two arbitrary regions of the real Universe. Phase randomization tests whether the marked statistic responds to mode relationships beyond the fixed two-point amplitudes.

In the sign-paired test, evaluate each statistic on \(\delta\) and \(-\delta\) and average. This removes all odd-order terms, including the third-order line above. The second-order line and the fourth-order normalization term cancel between matched-power fields. The disconnected Wick contractions within the four-field expectations also cancel because they are functions only of the matched two-point correlations.

For example, define the connected four-field projection \(T_{\delta s,\delta s}\) by

\[
\langle(\delta s)_1(\delta s)_2\rangle-u^2
=\xi_{\delta\delta}(r)\xi_{ss}(r)
+\xi_{\delta s}(r)\xi_{s\delta}(r)
+T_{\delta s,\delta s}(r).
\]

Thus the difference of the left-hand side equals \(\Delta T_{\delta s,\delta s}\) when the two-point functions match. The same Wick subtraction applies to the other four-field terms. For the \(p=2\) mark, the complete leading even-order difference is

\[
\boxed{\begin{aligned}
\Delta M(r)={}&
4\Delta T_{\delta s,\delta s}
+4\Delta T_{\delta s,s^2}
+\Delta T_{s^2,s^2}\\
&+2\Delta T_{\delta,\delta s^2}
+4\Delta T_{s,\delta s^2}
+O(\delta^6).
\end{aligned}}
\]

Here each \(\Delta T\) is the **injected-minus-control** connected four-point projection, not an assumption that the control's connected projection is exactly zero. The remainder starts at sixth order after sign pairing. At finite injection amplitude, it need not be negligible.

The code's centered quadratic operator \(q=\delta s-u\) isolates the first projection:

\[
\begin{aligned}
C_q(r)&=\langle q_1q_2\rangle,\\
T_q(r)&=C_q(r)-\xi_{\delta\delta}(r)\xi_{ss}(r)
-\xi_{\delta s}(r)\xi_{s\delta}(r),\\
\Delta C_q(r)&=\Delta T_q(r)
\quad\text{for exactly matched two-point functions}.
\end{aligned}
\]

Its contribution to \(\Delta M\) is therefore \(4\Delta T_q\) at fourth order, with **no fitted coefficient**. The difference \(\Delta M-4\Delta T_q\) contains the other four connected projections above and potentially higher even orders. The current 3D implementation bins individual two-point correlations before multiplying them in the Wick term; a paper-quality estimate should form Wick products at each separation vector before radial binning. The near-zero injected-minus-control identity error by itself checks matched-term cancellation, not the absolute accuracy of each binned \(T_q\).

## Linear-mark comparison and the three roles of the fields

With \(p=1\), the linear environmental mark is \(m=1+s\). Its marked contrast has denominator \(1+u\), and its cross-correlation with the unmarked density is

\[
X_E(r)=\langle\delta_1\delta_{m,2}\rangle
=\frac{\xi_{\delta\delta}(r)+\xi_{\delta s}(r)
+\langle\delta_1\delta_2s_2\rangle}{1+u}.
\]

The raw \(X_E\) is nonzero even when the bispectrum vanishes because its first two terms are two-point correlations. Under sign pairing, the cubic term cancels. For exactly matched-power fields, the remaining terms and denominator agree, so \(\Delta X_E=0\).

Three field roles must be kept distinct:

| Field | Role |
| --- | --- |
| \(g\) | Gaussian null, equivalent to injection amplitude \(A=0\) |
| \(\delta_A=g+A(g^3-3\langle g^2\rangle g)\) | Injected sign-symmetric field with zero ensemble bispectrum and a connected four-point response |
| Phase-randomized \(\delta_A\) | Artificial control with exactly matched Fourier amplitudes and power |

The odd cubic injection is related to kurtosis and the trispectrum, but the phase-randomized control is not literally a second Universe and is not guaranteed to have zero connected four-point signal in one finite realization. Multiple independent injection seeds and multiple phase surrogates per seed are needed for an uncertainty assessment.

## Why idealized CMB lensing is the physical example

Draw a Gaussian primary temperature field \(T_{\rm primary}\) independently of a Gaussian lensing potential \(\phi\), with \(C_\ell^{T\phi}=0\). Lensing remaps the primary temperature:

\[
T_{\rm lensed}(\hat n)
=T_{\rm primary}(\hat n+\nabla\phi(\hat n))
\]

where the expression is schematic on the sphere. At fixed \(\phi\), the remapping is linear in the Gaussian primary. Across realizations of \(\phi\), it couples temperature modes and produces a connected temperature trispectrum. Changing the primary \(T\) to \(-T\) changes the lensed map to its negative; the ensemble's odd temperature moments, including its bispectrum, vanish in this idealized independent-\(T,\phi\) setting.

For the direct algebraic CMB analogue, set \(\delta_T=(T-\langle T\rangle)/\sigma_T\) in the current toy benchmark and smooth it to obtain \(s\). Here \(\delta_T\) is a unit-variance temperature contrast, not the physical contrast \((T-\bar T)/\bar T\). Apply the *same* mark and normalization as in the density case:

\[
m(s)=(1+s)^2,\qquad
F=(1+\delta_T)m(s),\qquad
\delta_m^{\rm CMB}=F/\langle F\rangle-1.
\]

The factor \(1+\delta_T\) is an algebraic analogue of \(1+\delta\); CMB temperature is not a number-density field. The earlier CMB polynomial mark \(1+0.5s+0.5(s^2-\langle s^2\rangle)\), applied to \(T\) without this mean normalization, defines a different exploratory descriptor and does not test this \(M\).

Normalization must be applied separately to the two sign partners. Let
\[
x=\delta_T+2s+\delta_Ts^2,\quad y=2\delta_Ts+s^2,\quad
F_\pm=(1\pm\delta_T)(1\pm s)^2=1+y\pm x,\quad
\mu_\pm=1+\langle y\rangle\pm\langle x\rangle.
\]
With \(X=x-\langle x\rangle\) and \(Y=y-\langle y\rangle\), the exact sign-averaged spectrum is
\[
C_\ell^{\rm mark}=\frac12(\mu_+^{-2}+\mu_-^{-2})(C_\ell^{XX}+C_\ell^{YY})
+(\mu_+^{-2}-\mu_-^{-2})C_\ell^{XY}.
\]
For \(q=\delta_Ts-\langle\delta_Ts\rangle\) and \(r=s^2-\langle s^2\rangle\), \(Y=2q+r\). Therefore the coefficient of \(C_\ell^{qq}\) is exactly \(K_q=2(\mu_+^{-2}+\mu_-^{-2})\). The plotted connected component is \(K_qT_q\), evaluated for each map and each matched phase control before their difference. At unit variance, \(K_q\) can differ substantially from the weak-field value four. The remainder contains the \(q\)-\(r\) cross term, other contractions, and normalization effects; a nonzero component alone does not establish a resolved trispectrum response.

For each lensed realization, phase-randomized harmonic controls retain its measured \(|a_{\ell m}|\) and hence its \(C_\ell\). An unlensed Gaussian primary supplies a separate no-lensing null. A lensed-minus-unlensed difference alone is not pure higher-order response because lensing also changes the two-point spectrum. Real observed CMB maps can contain ISW–lensing correlations, foregrounds and other bispectrum sources, so they do not automatically satisfy this clean \(B=0\) premise.

The controlled question is: **with the temperature bispectrum absent, does the nonlinear marked descriptor respond to the injected lensing trispectrum beyond a matched-\(C_\ell\) null?** Demonstrating that response would establish trispectrum sensitivity in clean CMB maps. It would not show that the descriptor reconstructs \(\phi\), is optimal, or adds information beyond a lensing quadratic estimator; that comparison is reserved for a later study.

## What the ensemble test must show

Run many independent Gaussian seeds. In 3D, compare \(A=0\) and fixed nonzero \(A\); for CMB, compare unlensed and independent-\(T,\phi\) lensed simulations. Use sign-paired estimators to cancel odd terms, and several exact-power phase surrogates per injected realization. The surrogates are correlated controls, not independent injected universes.

Measure the unpaired collapsed three-point projection to confirm its ensemble null, the signed marked difference, the direct connected \(T_q\) difference, and the normalization-corrected \(K_qT_q\) component. Plot seed-based uncertainty bands and a separate Gaussian/unlensed null distribution. A nonzero single-realization curve, or the algebraic identity \(\Delta C_q=\Delta T_q\), is not a detection significance. This paper's target is a validated response to connected four-point structure at fixed marks, not a parameter scan or an information comparison with CMB lensing reconstruction.
