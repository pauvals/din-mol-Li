#!/bin/bash
#SBATCH --job-name=dana-1n
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:59:00
#SBATCH --output=dana-1node-%j.out

set -euo pipefail

repo=${SLURM_SUBMIT_DIR:?Submit this job from the DANA repository root}
source "$repo/benchmarks/serafin/load_environment.sh"

export DANA_BENCH_NODES=1
export DANA_CORES_PER_NODE=64
export DANA_BENCH_TILE_X=${DANA_BENCH_TILE_X:-24}
export DANA_BENCH_TILE_Y=${DANA_BENCH_TILE_Y:-24}
export DANA_BENCH_STEPS=${DANA_BENCH_STEPS:-100}
export DANA_BENCH_REPEATS=${DANA_BENCH_REPEATS:-3}
export DANA_BENCH_CONFIGS=${DANA_BENCH_CONFIGS:-"1:1 1:64 2:32 4:16 8:8 16:4 32:2 64:1"}
export DANA_BENCH_DIR=${DANA_BENCH_DIR:-$repo/results/serafin-$SLURM_JOB_ID-1node}

"$repo/benchmarks/serafin/run_scaling_slurm.sh"
