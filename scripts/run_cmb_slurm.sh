#!/usr/bin/env bash
#SBATCH --job-name=cmb-marked
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=cmb-marked-%j.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

: "${COUNT:?Set COUNT to the requested number of paired skies}"
PYTHON_BIN="${PYTHON_BIN:-python}"
CONFIG="${CONFIG:-configs/cmb_marked_ensemble.json}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/cmb-matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"

args=(run --config "$CONFIG" --count "$COUNT")
if [[ -n "${OUTPUT_DIR:-}" ]]; then
  args+=(--output-dir "$OUTPUT_DIR")
fi
"$PYTHON_BIN" -m cmb.marked_ensemble_pilot "${args[@]}"
