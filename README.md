# Marked-correlation pilots

## GitHub and cluster run

The repository tracks source code, fixed configurations, and methodology
documents. Generated run directories are ignored by Git; copy or archive
results separately. The CMB environment used locally was Python 3.11 with
the versions in [requirements-cmb.txt](requirements-cmb.txt). On an older
Linux cluster, pip may try to compile those exact versions from source.
Use the binary Conda packages in [environment-hpc.yml](environment-hpc.yml):

    conda env create -f environment-hpc.yml
    conda activate marked-cmb
    python -m pip install --no-deps lenspyx==2.0.52
    python -c "import camb, healpy, lenspyx, numpy, scipy; print('ready')"

The exact-pinned requirements file records the local software versions; it
is not the recommended installer for a cluster with an older Linux runtime.

The ensemble configuration is in
[configs/cmb_marked_ensemble.json](configs/cmb_marked_ensemble.json). Its
default count is only a placeholder; choose the count when submitting.
For a Slurm cluster, submit from the repository root after adapting time,
memory, and account settings in [scripts/run_cmb_slurm.sh](scripts/run_cmb_slurm.sh):

    sbatch --export=ALL,COUNT=20,PYTHON_BIN="$CONDA_PREFIX/bin/python",OUTPUT_DIR="$SCRATCH/cmb_marked_ensemble" scripts/run_cmb_slurm.sh

If the job reaches its wall-time limit, resubmit with the same count and
output directory; completed pairs are skipped. A larger count extends the
same deterministic seed sequence. Run status or diagnose with the same
--config, --count, and --output-dir arguments. The runner needs no plotting
backend; diagnose generates the figure and vector file later.

For the CMB setup, start with [cmb_current_configuration.md](cmb_current_configuration.md). The marked-spectrum filter pilots are in [run_cmb_marked_pilot_2048](run_cmb_marked_pilot_2048) and [run_cmb_marked_pilot_lmin600_2048](run_cmb_marked_pilot_lmin600_2048). The earlier angular and harmonic consistency check is in [run_cmb_angular_consistency_2048](run_cmb_angular_consistency_2048). Older outputs are catalogued in [archive/cmb_pilots](archive/cmb_pilots/README.md).

Active code lives in `cmb/` and `lss/`. Run it from the repository root with `python -m`, for example `python -m cmb.analyze_cmb_output_binning run_cmb_angular_consistency_2048`. The historical CMB demo and analysis scripts are in `archive/cmb_pilots/scripts/`.

`lss/lognormal_marked_correlation.py` implements the toy experiment described in
`nonlinear_marked_correlations_as_joint_compressed_probes_of_the_bispectrum_and_trispectrum.md`.
It requires NumPy. Matplotlib is optional for the diagnostic figure.

```bash
python3 -m lss.lognormal_marked_correlation --nmesh 128 --boxsize 500 \
  --smoothing-radius 10 --mark-power 2 --output-dir results
```

The script first generates a Gaussian log-density field with a smooth toy
spectrum, then constructs

\[
1+\delta=\exp(g-\mathrm{Var}(g)/2).
\]

It creates a control by phase-randomizing the *final* lognormal density field.
That preserves the amplitude of every discrete Fourier mode, hence its power
spectrum and two-point correlation, while destroying the phase correlations
responsible for the field's higher-order spatial structure.

For both fields it applies Gaussian harmonic-space smoothing,

\[
\delta_R(\mathbf{k})=e^{-k^2R^2/2}\delta(\mathbf{k}),
\]

forms the normalized marked fluctuation from

\[
F=(1+\delta)(1+\delta_R)^p,
\]

and estimates the periodic marked auto-correlation
\(M(r)=1+\xi_{F/\langle F\rangle-1}(r)\).  The principal diagnostic is
\(\Delta M=M_{\rm lognormal}-M_{\rm phase-randomized}\).

Outputs in the selected directory:

- `power_spectra.csv`: binned power spectra, verifying the matched-power control;
- `marked_correlations.csv`: ordinary and marked correlations plus \(\Delta M\);
- `marked_correlation_diagnostics.png`: the corresponding quick-look figure,
  when Matplotlib is installed.

The `--sigma-g` parameter controls lognormal non-Gaussianity.  `--mark-power`
and `--smoothing-radius` are the mark parameters to scan.  Keep the latter
larger than a few mesh cells to avoid interpreting grid-scale artifacts.

## Trispectrum-only comparison

`lss/trispectrum_marked_comparison.py` is the complementary first comparison with
Ebina's density--linear-mark cross-statistic.  It uses an odd cubic transform
of a Gaussian field, then averages the measurements from \(\delta\) and
\(-\delta\).  This sign pairing cancels every odd-order contribution exactly:
the paired field has \(B=0\), while its connected trispectrum is non-zero.
The nonlinear marked auto-correlation can therefore retain a residual after
matching the power spectrum, whereas the Ebina cross residual is zero.

```bash
python3 -m lss.trispectrum_marked_comparison --nmesh 64 --output-dir run_trispectrum_64
```

The output directory contains `trispectrum_comparison.csv`, matched power
spectra, and `trispectrum_comparison.png`.  The middle panel overlays the
nonlinear and Ebina statistics as correlation excesses with a logarithmic
separation axis; the last panel shows the nonlinear-auto residual magnitude in
log--log form.

## CMB ensemble benchmark

### Ten-sky marked-spectrum pilot

The current validation configuration is in [configs/cmb_marked_10.json](configs/cmb_marked_10.json).
It uses ten new paired lensed/unlensed skies at Nside 2048, input temperature
multipoles 600–3000, a 2 arcmin Gaussian smoothing sigma, the sign-paired
p=2 marked spectrum, and three matched phase controls per map. Output bins
cover 40 ≤ L < 1200. The diagnostic tracks 100 ≤ L < 600 and
600 ≤ L < 1200 separately.

Run from the repository root:

    env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MPLCONFIGDIR=/private/tmp/mpl_cmb \
      /opt/anaconda3/envs/pylians/bin/python -m cmb.marked_ensemble_pilot run

The command prints progress after each pair and is safe to rerun: completed
pairs are skipped. In another terminal, inspect the current result with

    /opt/anaconda3/envs/pylians/bin/python -m cmb.marked_ensemble_pilot status

or regenerate the plot with action diagnose. Results appear under
run_cmb_marked_10_2048: one folder and run.log per paired sky, plus summary.csv,
progress.json, and marked_diagnostic.png. The plotted lines are individual
skies and their mean; no uncertainty bars are inferred from this pilot.
The previous three-sky runtime suggests roughly 45–60 minutes for ten skies,
subject to the wider output range and machine load. This is a screening run
before selecting a larger ensemble.

The descriptor code now limits the quadratic and marked-field harmonic
transforms to the requested output multipoles and reuses the Gaussian Wick
term across matched phase controls. The input temperature harmonics still
extend to the configured lmax. A repeat of master seed 340001 took 214 seconds
versus 338 seconds in the original run. Its paired marked bandpowers changed
by less than 5e-13 in both summary bands. The optimized run records setup
and per-stage timings in benchmark.json; this one spent 148 seconds on six
phase controls, 46 seconds on the two main descriptors, 13 seconds filtering,
and 6 seconds generating maps. Adding phase controls therefore remains the
main runtime cost.

### Extendable CMB ensemble

[configs/cmb_marked_ensemble.json](configs/cmb_marked_ensemble.json) fixes the
same input and 43 output bins for a future independent ensemble. It starts
with ten phase controls per lensed or unlensed map and a seed generator whose
first N seeds do not change when the requested count grows. This uses a new
output directory and independent seeds, separate from the completed pilot.
Choose the number of paired skies at run time:

    env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 \
      /opt/anaconda3/envs/pylians/bin/python -m cmb.marked_ensemble_pilot run \
      --config configs/cmb_marked_ensemble.json --count 20

The run writes each pair to its own directory and updates summary.csv and
progress.json after completion. Restart with the same command to skip
completed pairs. A later larger --count extends the same seed sequence.
The runner does not import Matplotlib or generate a plot. Use the same
--config and --count options with action status for progress, or action
diagnose to create marked_diagnostic.png and ensemble_vectors.npz.
The latter contains the complete 43-bin lensed, unlensed, paired marked,
connected-q, and remainder vectors plus their sample marked covariance.
With fewer than 44 independent pairs that full-bin covariance is singular;
retain the vectors but use a smaller prespecified band set for covariance
based inference. For the decomposition and broad-band report, run
python -m cmb.analyze_marked_ensemble on the output directory.

The measured optimized cost with three controls was 214 seconds per pair.
Ten controls should take longer because each added control requires another
marked-field calculation; benchmark a small number of pairs before projecting
the time for a large ensemble.

### Compute-node notebook

The runner accepts --workers to process independent paired skies in separate
processes. Each process uses the config's four harmonic-transform threads.
On a large CPU node, start with eight workers (32 requested threads) and
monitor memory before increasing that count. Each worker writes its own
pair directory and run.log; the parent updates progress.json atomically.
The output directory can live on a scratch filesystem.

In a notebook using the CMB environment's Python kernel, set the repository
and scratch paths and run:

    import os, sys, subprocess
    from pathlib import Path

    repo = Path.home() / "MarkedCorrelationNL_Trispectrum"
    scratch = Path(os.environ["SCRATCH"])
    output = scratch / "cmb_marked_ensemble"
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "4"
    subprocess.run([
        sys.executable, "-m", "cmb.marked_ensemble_pilot", "run",
        "--config", "configs/cmb_marked_ensemble.json",
        "--count", "100", "--workers", "8",
        "--output-dir", str(output),
    ], cwd=repo, env=env, check=True)

The same call resumes completed pairs. Use action status with the same count
and output directory to inspect progress; run action diagnose after completion.
The notebook cell remains busy until the runner completes. Its own output is
short; individual simulation logs are in output/pair_XX/run.log.

For diagnostics in a notebook using this environment, no subprocess is
needed. Add the repository to Python's module path if the notebook is in
another directory, then call the analysis functions:

    import sys
    sys.path.insert(0, str(repo))
    from cmb.marked_ensemble_pilot import load_config, diagnose
    from cmb.analyze_marked_ensemble import analyze
    config = load_config((repo / "configs/cmb_marked_ensemble.json").resolve(),
                         count=100, output_dir=output)
    records = diagnose(config)
    print(analyze(output))

The benchmark creates paired lensed and unlensed maps from the same Gaussian
primary, and two phase controls for each map. It computes sign-paired marked
and connected-four-point bandpowers using cached spherical-harmonic transforms;
full maps are discarded after recording the bandpowers. The input spectra are
toy spectra, suitable for timing and software validation only.

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 \
  /opt/anaconda3/envs/pylians/bin/python -m cmb.benchmark_cmb_ensemble \
  --nside 128 --lmax 128 --nsims 10 --surrogates 2 \
  --nbins 12 --output-dir run_cmb_benchmark_128
```

The output contains `bandpowers.csv` (one row per multipole bin, seed, and
lensed/unlensed map) and `benchmark.json` (settings, timing, bin edges,
and the largest realized phase-control power mismatch). The ten-pair local
run took about 2.3 seconds, or 0.23 seconds per independent pair, at this
resolution. Rebenchmark if the spectrum range, resolution, or number of
phase controls changes.
