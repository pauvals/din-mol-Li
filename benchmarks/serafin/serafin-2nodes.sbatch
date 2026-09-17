#!/bin/bash
#SBATCH --job-name=dana-2n
#SBATCH --partition=multi
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=dana-2nodes-%j.out

set -euo pipefail

repo=${SLURM_SUBMIT_DIR:?Submit this job from the DANA repository root}
source "$repo/benchmarks/serafin/load_environment.sh"

export DANA_BENCH_NODES=2
export DANA_CORES_PER_NODE=64
export DANA_BENCH_TILE_X=${DANA_BENCH_TILE_X:-24}
export DANA_BENCH_TILE_Y=${DANA_BENCH_TILE_Y:-24}
export DANA_BENCH_STEPS=${DANA_BENCH_STEPS:-100}
export DANA_BENCH_REPEATS=${DANA_BENCH_REPEATS:-3}
export DANA_BENCH_CONFIGS=${DANA_BENCH_CONFIGS:-"2:64 4:32 8:16 16:8 32:4 64:2 128:1"}
export DANA_BENCH_DIR=${DANA_BENCH_DIR:-$repo/results/serafin-$SLURM_JOB_ID-2nodes}

"$repo/benchmarks/serafin/run_scaling_slurm.sh"
