# Archived CMB pilots

These runs document exploratory choices made before the current CMB configuration. The current full-sky pilot with raw spectra and the angular/harmonic consistency check remains in run_cmb_angular_consistency_2048/.

| Run | Historical purpose |
| --- | --- |
| run_cmb_lensing_64/ and run_cmb_marked_64/ | Low-resolution lensing demo and a polynomial mark distinct from the current \(p=2\) mark |
| run_cmb_benchmark_128/ and run_cmb_benchmark_128_p2/ | Low-resolution benchmarks preceding the current input temperature and smoothing choices |
| run_cmb_actlike_single_2048/ and run_cmb_actlike_paired_controls_2048/ | First full-resolution pair and phase-control pairing check |
| run_cmb_actlike_fine_2048/ | Restricted \(L<200\) output with finer bins |
| run_cmb_actlike_high_ell_2048/ | Extended output for the input temperature cut \(\ell\geq600\) |
| run_cmb_actlike_lmin100_2048/ | Earlier copy of the \(\ell\geq100\) pair, superseded by the active angular-consistency run with saved raw spectra |

The benchmark.json files in moved runs now show their archive location as output_dir and retain the generation path as original_output_dir.

Historical analysis scripts and the earlier CMB exercise recap are also kept here. Run an archived script from the repository root with the repository on the Python import path when it imports a current root module.
